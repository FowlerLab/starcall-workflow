'use strict';

// ---------- state ----------
const state = {
  config: null,
  mode: 'dot',
  data: null,          // current /api/dot or /api/cell response
  images: [],          // every thumbnail spec on the page, in prev/next order
  selected: -1,
  overrides: {},       // 'seq:G' -> {vmin, vmax}, set with "Apply to channel"
  viewerRange: null,   // unsaved slider range of the image in the viewer
};
const $ = (id) => document.getElementById(id);
const el = (tag, attrs = {}, ...children) => {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === 'class') node.className = v;
    else if (k === 'text') node.textContent = v;
    else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v);
  }
  for (const child of children) if (child != null) node.append(child);
  return node;
};
const SVG_NS = 'http://www.w3.org/2000/svg';
const svg = (tag, attrs = {}) => {
  const node = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  return node;
};
const fmt = (v, digits = 3) => (v == null || Number.isNaN(v)) ? '–' : Number(v).toFixed(digits);

// ---------- channel colours ----------
// Each sequencing channel has one colour, used for its image column, the chart lines and the
// called-base chips. Bases default to the chart series slots; DAPI and GFP to blue and green.
// Phenotype channels are coloured separately, keyed 'pt:<channel>'.
const BASE_ORDER = ['G', 'T', 'A', 'C'];
const FIXED_DEFAULTS = { DAPI: '#3d5afe', GFP: '#22c32e' };
const PT_DEFAULTS = { DAPI: '#3d5afe', GFP: '#22c32e', 'Ph+WGA': '#ff4d4d', Mito: '#e040fb' };
const colorKey = (image) => image.pt ? 'pt:' + image.channel : image.channel;
const COLOR_KEY = 'starcall-viewer-colors';
state.colors = (() => {
  try { return JSON.parse(localStorage.getItem(COLOR_KEY)) || {}; } catch (e) { return {}; }
})();
function saveColors() {
  try { localStorage.setItem(COLOR_KEY, JSON.stringify(state.colors)); } catch (e) { /* storage unavailable */ }
}
function defaultColor(channel) {
  if (channel.startsWith('pt:')) return PT_DEFAULTS[channel.slice(3)] || '#ffffff';
  if (channel in FIXED_DEFAULTS) return FIXED_DEFAULTS[channel];
  const i = BASE_ORDER.indexOf(channel);
  if (i >= 0) return getComputedStyle(document.documentElement).getPropertyValue(`--series-${i + 1}`).trim();
  return '#ffffff';
}
const channelColor = (channel) => state.colors[channel] || defaultColor(channel);

function colorPicker(key, name) {
  const picker = el('input', { type: 'color', value: channelColor(key), 'data-color-channel': key,
    title: `Colour for ${name} images`, 'aria-label': `Colour for ${name}` });
  picker.addEventListener('input', () => { state.colors[key] = picker.value; saveColors(); applyColors(); });
  return picker;
}
const baseColor = (base) => BASE_ORDER.includes(base) ? channelColor(base) : 'var(--axis)';

function hexColor(color) {
  // '#rgb' / '#rrggbb' -> 'rrggbb' for the server's colour map
  let hex = String(color).trim().replace(/^#/, '');
  if (/^[0-9a-f]{3}$/i.test(hex)) hex = [...hex].map((c) => c + c).join('');
  return /^[0-9a-f]{6}$/i.test(hex) ? hex.toLowerCase() : 'ffffff';
}

function refreshImages() {
  // re-request every image already on the page (new colour or contrast), keeping the same elements
  for (const image of state.images) {
    const url = cropUrl(image);
    if (image.img.getAttribute('src') !== url) image.img.src = url;
  }
  const image = state.images[state.selected];
  if (image && state.viewerImg) {
    const url = cropUrl(image, state.viewerRange);
    if (state.viewerImg.getAttribute('src') !== url) state.viewerImg.src = url;
  }
}

function applyColors() {
  // images are re-rendered by the server in the new colour; wait until the picker settles
  clearTimeout(applyColors.timer);
  applyColors.timer = setTimeout(refreshImages, 150);
  document.querySelectorAll('[data-base]').forEach((node) => { node.style.borderBottomColor = baseColor(node.dataset.base); });
  document.querySelectorAll('input[data-color-channel]').forEach((input) => { input.value = channelColor(input.dataset.colorChannel); });
  if (state.data && state.data.kind === 'dot') renderDotCharts();
}

// ---------- api ----------
async function api(path) {
  const response = await fetch(path);
  if (!response.ok) {
    let body = {};
    try { body = await response.json(); } catch (e) { /* not json */ }
    const error = new Error(typeof body.detail === 'string' ? body.detail : response.statusText);
    if (body.expected) { error.expected = body.expected; error.missing = body.missing; }
    throw error;
  }
  return response.json();
}

function contrastFor(image) {
  const key = `${image.pt ? 'pt' : 'seq'}:${image.channel}`;
  return state.overrides[key] || null;
}

function cropUrl(image, range) {
  if (image.composite) return compositeUrl();
  const d = state.data;
  const box = image.pt ? d.box_pt : d.box;
  const params = new URLSearchParams({ well: state.well, tile: d.tile, cycle: image.cycle, channel: image.channel, box: box.join(',') });
  range = range || contrastFor(image);
  if (range) { params.set('vmin', range.vmin); params.set('vmax', range.vmax); }
  else params.set('contrast', $('contrast-mode').value);
  params.set('color', hexColor(channelColor(colorKey(image))));
  return '/api/crop.png?' + params;
}

// ---------- composite ----------
// Images stacked into one, each from black to its channel colour with its contrast, added together.
// Layers are keyed 'pt|PT|DAPI' / 'seq|01|G'; the choice is remembered in the browser, and until
// one is made the composite is the phenotype segmentation channels of the PT cycle.
const COMPOSITE_KEY = 'starcall-viewer-composite';
const EMPTY_IMAGE = 'data:image/gif;base64,R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7';
const layerKey = (layer) => `${layer.pt ? 'pt' : 'seq'}|${layer.cycle}|${layer.channel}`;
state.composite = (() => {
  try { return JSON.parse(localStorage.getItem(COMPOSITE_KEY)); } catch (e) { return null; }
})();
function saveComposite() {
  try {
    if (state.composite) localStorage.setItem(COMPOSITE_KEY, JSON.stringify(state.composite));
    else localStorage.removeItem(COMPOSITE_KEY);
  } catch (e) { /* storage unavailable */ }
}

function defaultLayerKeys() {
  const c = state.config;
  const cycle = c.pt_cycles.includes('PT') ? 'PT' : c.pt_cycles[0];
  return cycle ? c.segmentation_channels.map((channel) => layerKey({ pt: true, cycle, channel })) : [];
}
const selectedLayerKeys = () => new Set(Array.isArray(state.composite) ? state.composite : defaultLayerKeys());

function availableLayers() {
  // every image the well has, phenotype first, in the order of the page
  const c = state.config;
  const layers = [];
  for (const cycle of c.pt_cycles) for (const channel of c.phenotyping_channels) layers.push({ pt: true, cycle, channel });
  for (const cycle of c.seq_cycles) for (const channel of c.sequencing_channels) layers.push({ pt: false, cycle, channel });
  return layers;
}
function compositeLayers() {
  const keys = selectedLayerKeys();
  return availableLayers().filter((layer) => keys.has(layerKey(layer)));
}
const layerName = (layer) => layer.pt ? `${layer.cycle} ${layer.channel}` : `cycle ${layer.cycle} ${layer.channel}`;

function compositeLabel() {
  const layers = compositeLayers();
  if (!layers.length) return 'Composite (no layers)';
  if (layers.length > 4) return `Composite of ${layers.length} images`;
  return 'Composite · ' + layers.map(layerName).join(' + ');
}

function compositeUrl() {
  const layers = compositeLayers();
  if (!layers.length) return EMPTY_IMAGE;
  const d = state.data;
  const params = new URLSearchParams({ well: state.well, tile: d.tile, box: d.box.join(','), contrast: $('contrast-mode').value });
  for (const layer of layers) {
    const parts = [layer.cycle, layer.channel, hexColor(channelColor(colorKey(layer)))];
    const range = contrastFor(layer);
    if (range) parts.push(range.vmin, range.vmax);
    params.append('layer', parts.join('|'));
  }
  return '/api/composite.png?' + params;
}

function renderComposite() {
  // the composite image, after every other image so it doesn't change their prev/next order
  const image = { composite: true, label: compositeLabel() };
  $('composite-stage').replaceChildren(makeThumb(image));
  renderLayerPicker();
}

function renderLayerPicker() {
  const c = state.config;
  const keys = selectedLayerKeys();
  const group = (title, cycles, channels, pt) => {
    if (!cycles.length) return null;
    const rows = cycles.flatMap((cycle) => channels.map((channel) => ({ pt, cycle, channel })));
    const count = rows.filter((layer) => keys.has(layerKey(layer))).length;
    const table = el('table', { class: 'layer-table' });
    table.append(el('thead', {}, el('tr', {}, el('th', { text: pt ? '' : 'cycle' }),
      ...channels.map((channel) => el('th', {}, colorPicker(pt ? 'pt:' + channel : channel, pt ? `phenotype ${channel}` : channel), channel)))));
    const body = el('tbody');
    for (const cycle of cycles) {
      const row = el('tr', {}, el('th', { text: cycle }));
      for (const channel of channels) {
        const layer = { pt, cycle, channel };
        const box = el('input', { type: 'checkbox', 'aria-label': layerName(layer), title: layerName(layer) });
        box.checked = keys.has(layerKey(layer));
        box.addEventListener('change', () => toggleLayer(layerKey(layer), box.checked));
        row.append(el('td', {}, box));
      }
      body.append(row);
    }
    table.append(body);
    const details = el('details', { class: 'layer-group' },
      el('summary', {}, title, ' ', el('span', { class: 'count', text: `(${count} of ${rows.length})` })),
      el('div', { class: 'layer-table-wrap' }, table));
    // phenotype open; sequencing (many cycles) only when one of its images is in the composite
    details.open = pt || count > 0;
    return details;
  };
  $('composite-layers').replaceChildren(...[
    group('Phenotype', c.pt_cycles, c.phenotyping_channels, true),
    group('Sequencing', c.seq_cycles, c.sequencing_channels, false),
  ].filter(Boolean));
}

function toggleLayer(key, on) {
  const keys = selectedLayerKeys();
  if (on) keys.add(key); else keys.delete(key);
  state.composite = [...keys];
  saveComposite();
  updateComposite();
  for (const details of $('composite-layers').querySelectorAll('details')) {
    const boxes = details.querySelectorAll('input[type=checkbox]');
    details.querySelector('.count').textContent = `(${[...boxes].filter((b) => b.checked).length} of ${boxes.length})`;
  }
}

function updateComposite() {
  const image = state.images.find((other) => other.composite);
  if (!image) return;
  image.label = compositeLabel();
  image.node.title = image.img.alt = image.label;
  const empty = !compositeLayers().length;
  document.querySelectorAll('.composite-empty').forEach((tag) => { tag.hidden = !empty; });
  const url = cropUrl(image);
  if (image.img.getAttribute('src') !== url) { image.img.classList.add('loading-img'); image.img.src = url; }
  if (state.images[state.selected] === image && state.viewerImg) {
    $('viewer-title').textContent = image.label;
    state.viewerImg.alt = image.label;
    if (state.viewerImg.getAttribute('src') !== url) { state.viewerImg.classList.add('loading-img'); state.viewerImg.src = url; }
  }
}

function pathLayer(outline, className, stroke, width) {
  // an outline path along mask pixel edges (from the server), drawn at a fixed screen width
  if (!outline || !outline.path) return null;
  const [rows, cols] = outline.size;
  const layer = svg('svg', { class: 'overlay ' + className, viewBox: `0 0 ${cols} ${rows}`, preserveAspectRatio: 'none' });
  layer.append(svg('path', { d: outline.path, fill: 'none', stroke, 'stroke-width': width,
    'stroke-linecap': 'square', 'vector-effect': 'non-scaling-stroke' }));
  layer.style.pointerEvents = 'none';
  return layer;
}

function outlineLayer() {
  // the view's cell
  return pathLayer(state.data.cell && state.data.cell.outline, 'outline', 'var(--outline)', 1.5);
}

function nucleusLayer() {
  // the nucleus matched to the view's cell
  return pathLayer(state.data.cell && state.data.cell.nucleus_outline, 'nucleus', 'var(--nucleus)', 1.5);
}

// ---------- routing ----------
function parseHash() {
  const [path, query] = decodeURIComponent(location.hash.replace(/^#\/?/, '')).split('?');
  const parts = path.split('/').filter(Boolean);
  const mode = parts[0] === 'cell' ? 'cell' : 'dot';
  // links from before wells were selectable have no well part
  const well = parts[1] && state.config.wells.includes(parts[1]) ? parts.splice(1, 1)[0] : null;
  if (mode === 'dot' && parts.length >= 3) {
    const size = new URLSearchParams(query || '').get('size');
    return { mode, well, tile: parts[1], index: parseInt(parts[2], 10), size: size ? parseInt(size, 10) : null };
  }
  if (mode === 'cell' && parts.length >= 2) {
    const size = new URLSearchParams(query || '').get('size');
    return { mode, well, id: parseInt(parts[1], 10), size: size ? parseInt(size, 10) : null };
  }
  return { mode, well };
}

const dotHash = (tile, index, size) => `#dot/${state.well}/${tile}/${index}` + (size && size !== 'fit' ? `?size=${size}` : '');
const cellHash = (id, size) => `#cell/${state.well}/${id}` + (size && size !== 'fit' ? `?size=${size}` : '');

function setMode(mode) {
  state.mode = mode;
  $('tab-dot').setAttribute('aria-selected', mode === 'dot');
  $('tab-cell').setAttribute('aria-selected', mode === 'cell');
  $('dot-form').hidden = mode !== 'dot';
  $('cell-form').hidden = mode !== 'cell';
}

async function route() {
  const target = parseHash();
  setMode(target.mode);
  const well = target.well || state.well || state.config.default_well;
  if (well !== state.well || !state.wellReady) {
    if (!(await selectWell(well))) return;
  }
  if (target.mode === 'dot' && target.tile != null && !Number.isNaN(target.index)) {
    // 'whole' (an untiled run) is used as is; short grid names like '02x03y' get the 'tile' prefix
    const tile = state.config.tiles.includes(target.tile) || target.tile.startsWith('tile') ? target.tile : 'tile' + target.tile;
    $('dot-tile').value = tile;
    $('dot-index').value = target.index;
    $('dot-size').value = target.size || 'fit';
    if (!$('dot-size').value) $('dot-size').value = 'fit';   // a size from a link that isn't in the list
    const params = new URLSearchParams({ well: state.well, tile, index: target.index });
    if (target.size) params.set('size', target.size);
    await load(`/api/dot?${params}`);
  } else if (target.mode === 'cell' && !Number.isNaN(target.id) && target.id != null) {
    $('cell-id').value = target.id;
    $('cell-size').value = target.size || 'fit';
    const params = new URLSearchParams({ well: state.well, id: target.id });
    if (target.size) params.set('size', target.size);
    await load(`/api/cell?${params}`);
  } else {
    clearView();
  }
}

async function selectWell(well) {
  // per-well cycles, tiles and contrast; reports missing files if the well can't be shown
  state.well = well;
  state.wellReady = false;
  $('well').value = well;
  $('error').hidden = true;
  try {
    const info = await api(`/api/well?${new URLSearchParams({ well })}`);
    Object.assign(state.config, info);
  } catch (error) {
    clearView();
    showError(`Can't show ${well}`, error);
    return false;
  }
  for (const id of ['dot-tile', 'cell-tile']) {
    const select = $(id), current = select.value;
    select.replaceChildren(...state.config.tiles.map((tile) => el('option', { value: tile, text: tile })));
    if (state.config.tiles.includes(current)) select.value = current;
    // an untiled run has one tile, the whole well, so there is nothing to choose
    select.closest('label').hidden = !state.config.tiled;
  }
  state.wellReady = true;
  return true;
}

function clearView() {
  state.data = null;
  state.images = [];
  $('view').hidden = true;
  $('placeholder').hidden = false;
}

function showError(title, error) {
  const box = $('error');
  box.replaceChildren(el('div', { class: 'error-title', text: `${title}: ${error.message}` }));
  if (error.expected) {
    const fileList = (files) => el('ul', { class: 'file-list' }, ...files.map((f) =>
      el('li', { class: f.exists ? 'present' : 'absent' },
        el('span', { class: 'mark', text: f.exists ? '✓' : '✗', 'aria-label': f.exists ? 'present' : 'missing' }),
        el('span', { text: f.label }), el('code', { text: f.path }))));
    box.append(
      el('h3', { text: `Missing (${error.missing.length})` }), fileList(error.missing),
      el('details', {}, el('summary', { text: `All expected files (${error.expected.length})` }), fileList(error.expected)));
  }
  box.hidden = false;
  $('placeholder').hidden = true;
}

// ---------- loading ----------
async function load(url) {
  $('placeholder').hidden = true;
  $('error').hidden = true;
  $('view').hidden = false;
  $('view').className = 'view-' + state.mode;
  $('charts').hidden = state.mode !== 'dot';
  $('summary-images').hidden = true;
  $('summary').replaceChildren(el('div', { class: 'loading', text: 'Loading… (the first view of a tile builds its table cache, which can take a minute)' }));
  for (const id of ['pt-grid', 'matrix', 'charts', 'viewer-stage', 'cell-extra', 'composite-stage', 'composite-layers']) $(id).replaceChildren();
  $('viewer-title').textContent = 'Click any image to view it here';
  $('tooltip').hidden = true;
  try {
    const data = await api(url);
    state.data = data;
    state.region = null;
    // the cell id decides the cell; the tile dropdown follows it
    if (data.kind === 'cell') $('cell-tile').value = data.tile;
    render();
  } catch (error) {
    $('view').hidden = true;
    showError('Could not load', error);
  }
}

// ---------- rendering ----------
function render() {
  const d = state.data, c = state.config;
  state.images = [];
  $('view').className = 'view-' + d.kind;
  renderSummary();

  // cell view: segmentation channels in the top card
  const summaryImages = $('summary-images');
  summaryImages.replaceChildren();
  summaryImages.hidden = d.kind !== 'cell';
  if (d.kind === 'cell') {
    for (const cycle of c.pt_cycles) {
      for (const channel of c.segmentation_channels) {
        const image = { cycle, channel, pt: true, label: `${cycle} · ${channel} (segmentation)` };
        summaryImages.append(el('figure', {}, makeThumb(image),
          el('figcaption', {}, colorPicker(colorKey(image), `phenotype ${channel}`), `${cycle} · ${channel}`)));
      }
    }
  }

  // phenotype panel: segmentation channels for dots, every phenotype channel for cells
  const ptChannels = d.kind === 'dot' ? c.segmentation_channels : c.phenotyping_channels;
  $('pt-title').textContent = d.kind === 'dot' ? 'Phenotype (segmentation channels only)' : 'Phenotype (all channels)';
  const ptGrid = $('pt-grid');
  for (const cycle of c.pt_cycles) {
    for (const channel of ptChannels) {
      const image = { cycle, channel, pt: true, label: `${cycle} · ${channel}` };
      ptGrid.append(el('figure', {}, makeThumb(image),
        el('figcaption', {}, colorPicker(colorKey(image), `phenotype ${channel}`), image.label)));
    }
  }

  // sequencing matrix: rows are cycles, columns are config sequencing_channels in order
  const matrix = $('matrix');
  const channels = c.sequencing_channels;
  // the row label column is sized to its text by fitRowHeads(), and mirrored as padding on the right
  matrix.style.gridTemplateColumns = `var(--row-head) repeat(${channels.length}, minmax(64px, 1fr))`;
  matrix.append(el('div', { class: 'corner-head', text: 'cycle' }));
  for (const channel of channels) {
    matrix.append(el('div', { class: 'col-head' }, colorPicker(channel, channel), el('span', { text: channel })));
  }
  c.seq_cycles.forEach((cycle, i) => {
    const head = el('div', { class: 'row-head', title: `cycle ${cycle}` }, cycle);
    if (d.kind === 'dot' && d.sequence && d.sequence[i]) {
      const b = el('b', { text: d.sequence[i], title: 'called base', 'data-base': d.sequence[i] });
      b.style.borderBottomColor = baseColor(d.sequence[i]);
      head.append(b);
    }
    matrix.append(head);
    for (const channel of channels) matrix.append(makeThumb({ cycle, channel, pt: false, label: `cycle ${cycle} · ${channel}` }));
  });

  fitRowHeads();
  renderComposite();

  if (d.kind === 'dot') renderDotCharts();
  else renderCellExtra();
  $('charts').hidden = d.kind !== 'dot';
  $('cell-extra').hidden = d.kind !== 'cell';

  // start with the first phenotype image in the viewer
  select(0);
  if ($('show-region').checked) loadRegion();
}

async function loadRegion() {
  const d = state.data;
  if (!d) return;
  const params = new URLSearchParams({ well: state.well, tile: d.tile, box: d.box.join(',') });
  if (d.cell) params.set('cell', d.cell.label);
  $('region-status').textContent = 'loading…';
  try {
    const region = await api('/api/region?' + params);
    if (state.data !== d || !$('show-region').checked) return;   // the view or the option changed meanwhile
    state.region = region;
    $('region-status').textContent = `${region.num_cells} other cells, ${region.dots.length} dots`;
    refreshOverlays();
  } catch (error) {
    if (state.data === d) $('region-status').textContent = 'failed: ' + error.message;
  }
}

function refreshOverlays() {
  for (const image of state.images) setOverlays(image.node, false);
  const stage = $('viewer-stage').firstElementChild;
  if (stage) setOverlays(stage, true);
}

function fitRowHeads() {
  // narrowest row label column that fits every label, so the matrix card has no spare padding
  const matrix = $('matrix');
  const heads = [...matrix.querySelectorAll('.row-head, .corner-head')];
  if (!heads.length) return;
  const range = document.createRange();
  const width = Math.max(...heads.map((head) => { range.selectNodeContents(head); return range.getBoundingClientRect().width; }));
  matrix.style.setProperty('--row-head', Math.ceil(width + 4) + 'px');
}

function renderSummary() {
  const d = state.data;
  const ids = el('div', { class: 'ids' });
  const idBlock = (label, value) => el('div', {}, el('div', { class: 'id-label', text: label }), el('div', { class: 'id-value', text: value }));
  const facts = el('div', { class: 'facts' });
  const summary = $('summary');

  if (d.kind === 'dot') {
    if (state.config.tiled) ids.append(idBlock('Tile', d.tile));
    ids.append(idBlock('Dot index', d.index));
    const seq = el('div', { class: 'sequence', 'aria-label': 'sequence ' + d.sequence });
    [...d.sequence].forEach((base, i) => {
      const chip = el('span', { title: `cycle ${d.cycles[i]}`, 'data-base': base }, base, el('small', { text: i + 1 }));
      chip.style.borderBottomColor = baseColor(base);
      seq.append(chip);
    });
    facts.append(el('span', { text: `mean chastity ${fmt(d.mean_chastity)}` }), el('span', { text: `min chastity ${fmt(d.min_chastity)}` }));
    if (d.cell) {
      facts.append(el('span', {}, 'cell ', el('a', { href: cellHash(d.cell.id), text: d.cell.id }), ` (label ${d.cell.label})`));
    } else {
      facts.append(el('span', { text: `not in a cell` }));
    }
    if (!d.cell || $('dot-size').value !== 'fit') {
      facts.append(el('span', { text: `showing ${d.box[2] - d.box[0]} px around the dot` }));
    }
    facts.append(...coverageNotes());
    summary.replaceChildren(ids, seq, facts);
  } else {
    ids.append(idBlock('Cell id', d.cell.id));
    if (state.config.tiled) ids.append(idBlock('Tile', d.tile));
    facts.append(
      el('span', { text: `label ${d.cell.label}` }),
      el('span', { text: `${d.num_reads ? parseFloat(d.num_reads) : 0} unique reads` }),
      el('span', { text: `${d.total_count ? parseFloat(d.total_count) : 0} dots in reads table` }),
      el('span', { text: `${d.dots.length} dots shown` }));
    facts.append(...coverageNotes());
    summary.replaceChildren(ids, facts);
  }
}

function makeThumb(image) {
  const d = state.data;
  const index = state.images.length;
  image.index = index;
  state.images.push(image);

  const node = el('div', { class: 'thumb', role: 'button', tabindex: '0', title: image.label });
  const img = el('img', { alt: image.label, loading: 'lazy', class: 'crop loading-img' });
  img.addEventListener('load', () => img.classList.remove('loading-img'));
  img.src = cropUrl(image);
  node.append(...[img, noDataTag(image)].filter(Boolean));
  setOverlays(node, false);
  node.addEventListener('click', () => select(index));
  node.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); select(index); } });
  image.node = node;
  image.img = img;
  return node;
}

function noDataTag(image) {
  if (image.composite) {
    const tag = el('span', { class: 'no-data composite-empty', text: 'no images ticked' });
    tag.hidden = compositeLayers().length > 0;
    return tag;
  }
  // images with no frame at all in the crop: the stage never imaged this spot in that cycle
  const coverage = (state.data.coverage || {})[image.cycle];
  if (coverage !== 0) return null;
  return el('span', { class: 'no-data', text: image.pt ? `not imaged in ${image.cycle}` : `not imaged in cycle ${image.cycle}` });
}

function coverageNotes() {
  // header notes on cycles whose images don't cover the whole crop
  const coverage = state.data.coverage || {};
  const cycles = [...state.config.seq_cycles, ...state.config.pt_cycles];
  const none = cycles.filter((c) => coverage[c] === 0);
  const partly = cycles.filter((c) => coverage[c] > 0 && coverage[c] < 0.999);
  const notes = [];
  if (none.length) notes.push(el('span', { class: 'coverage-note', text: `no image data in cycles ${none.join(', ')}` }));
  if (partly.length) notes.push(el('span', { class: 'coverage-note', text: `partly imaged in cycles ${partly.join(', ')}` }));
  return notes;
}

function setOverlays(node, interactive) {
  // (re)draws outlines, dot circles and, in the enlarged view, dot labels on top of an image
  node.querySelectorAll(':scope > .overlay').forEach((layer) => layer.remove());
  // every layer lets clicks through except its own click targets, so drawing order doesn't affect clicking
  const layers = [regionOutlineLayer(interactive), regionNucleusLayer(), outlineLayer(), nucleusLayer(), markerLayer(interactive)];
  if (interactive) layers.push(labelLayer());
  node.append(...layers.filter(Boolean));
}

function shownDots() {
  // the view's own dots, plus every other dot in the crop when 'All cells and dots in view' is on
  const d = state.data;
  if (!state.region) return d.dots;
  const own = new Set(d.dots.map((dot) => dot.index));
  return d.dots.concat(state.region.dots.filter((dot) => !own.has(dot.index)).map((dot) => ({ ...dot, region: true })));
}

function regionOutlineLayer(interactive = false) {
  // outlines of the other cells in the crop, under the view's own cell outline. In the enlarged
  // view each cell is its own path, and clicking one opens that cell
  if (!interactive || !state.region || !state.region.cells) {
    return pathLayer(state.region && state.region.outline, 'outline region', 'var(--outline-other)', 1);
  }
  const [rows, cols] = state.region.outline.size;
  const layer = svg('svg', { class: 'overlay outline region', viewBox: `0 0 ${cols} ${rows}`, preserveAspectRatio: 'none' });
  layer.style.pointerEvents = 'none';
  for (const cell of state.region.cells) {
    const group = svg('g', { class: 'region-cell' });
    const line = svg('path', { class: 'region-line', d: cell.outline.path, fill: 'none', stroke: 'var(--outline-other)',
      'stroke-width': 1, 'stroke-linecap': 'square', 'vector-effect': 'non-scaling-stroke' });
    // the cell's whole area takes the click: outlines of touching cells share edges, so a click
    // on a line alone could belong to either cell
    const hit = svg('path', { class: 'region-hit', d: cell.area, fill: 'transparent', stroke: 'none' });
    hit.style.pointerEvents = 'fill';
    hit.style.cursor = 'pointer';
    const title = svg('title');
    title.textContent = `cell ${cell.id} (click to open)`;
    hit.append(title);
    hit.addEventListener('click', () => { location.hash = cellHash(cell.id, $('cell-size').value); });
    group.append(hit, line);
    layer.append(group);
  }
  return layer;
}

function regionNucleusLayer() {
  // nuclei of the other cells in the crop
  return pathLayer(state.region && state.region.nuclei_outline, 'nucleus region', 'var(--nucleus-other)', 1);
}

function markerLayer(interactive = false) {
  const d = state.data;
  const side = d.box[2] - d.box[0];
  const layer = svg('svg', { class: 'overlay markers', viewBox: '0 0 1 1', preserveAspectRatio: 'none' });
  for (const dot of shownDots()) {
    const r = (dot.current ? 5 : 2.5) / side;
    const circle = svg('circle', {
      cx: dot.rel[1], cy: dot.rel[0], r, fill: 'none',
      stroke: dot.current ? 'var(--dot-ring)' : dot.region ? 'var(--dot-ring-region)' : 'var(--dot-ring-other)',
      'stroke-width': dot.current ? 2 : 1.5, 'vector-effect': 'non-scaling-stroke',
    });
    if (interactive && !dot.current) {
      circle.style.pointerEvents = 'all';
      circle.style.cursor = 'pointer';
      const title = svg('title');
      title.textContent = `dot ${dot.index} · ${dot.sequence} · min chastity ${fmt(dot.min_chastity)} (click to open)`;
      circle.append(title);
      circle.addEventListener('click', () => { location.hash = dotHash(d.tile, dot.index); });
    }
    layer.append(circle);
  }
  // only the circles themselves take clicks, so outlines underneath stay clickable
  layer.style.pointerEvents = 'none';
  return layer;
}

function labelLayer() {
  // dot ids and sequences next to each dot circle in the enlarged view, shown by the
  // 'Cells and dot IDs / sequences' modes (too crowded to read on the thumbnails)
  const d = state.data;
  const side = d.box[2] - d.box[0];
  const layer = el('div', { class: 'overlay labels' });
  for (const dot of shownDots()) {
    const r = (dot.current ? 5 : 2.5) / side;
    for (const [kind, text] of [['id', dot.index], ['seq', dot.sequence]]) {
      const label = el('span', { class: 'dot-label ' + kind, text });
      label.style.left = `calc(${(dot.rel[1] + r) * 100}% + 2px)`;
      label.style.top = `${dot.rel[0] * 100}%`;
      layer.append(label);
    }
  }
  return layer;
}

// ---------- viewer ----------
function select(index) {
  if (!state.images.length) return;
  index = (index + state.images.length) % state.images.length;
  const prev = state.images[state.selected];
  if (prev && prev.node) prev.node.classList.remove('selected');
  state.selected = index;
  const image = state.images[index];
  image.node.classList.add('selected');
  $('viewer-title').textContent = image.label;

  const stage = $('viewer-stage');
  const node = el('div', { class: 'thumb' });
  const img = el('img', { alt: image.label, class: 'crop loading-img' });
  img.addEventListener('load', () => img.classList.remove('loading-img'));
  img.src = cropUrl(image);
  node.append(...[img, noDataTag(image)].filter(Boolean));
  setOverlays(node, true);
  stage.replaceChildren(node);
  state.viewerImg = img;
  state.viewerRange = null;

  // a composite has no single range; each layer uses its own channel's
  for (const id of ['vmin', 'vmax', 'apply-contrast', 'reset-contrast']) $(id).disabled = !!image.composite;
  if (image.composite) {
    $('vmin-out').value = $('vmax-out').value = '';
    return;
  }
  const kind = image.pt ? 'pt' : 'seq';
  const defaults = (state.config.contrast[kind] || {})[image.channel] || { vmin: 0, vmax: 65535, max: 65535 };
  const range = contrastFor(image) || { vmin: defaults.vmin, vmax: defaults.vmax };
  const top = Math.max(defaults.max, range.vmax) * 1.05;
  for (const id of ['vmin', 'vmax']) { $(id).max = Math.ceil(top); $(id).value = range[id]; }
  updateRangeOutputs();
}

function updateRangeOutputs() {
  $('vmin-out').value = Math.round($('vmin').value);
  $('vmax-out').value = Math.round($('vmax').value);
}

function onSlider() {
  let vmin = parseFloat($('vmin').value), vmax = parseFloat($('vmax').value);
  if (vmax <= vmin) { vmax = vmin + 1; $('vmax').value = vmax; }
  updateRangeOutputs();
  state.viewerRange = { vmin, vmax };
  clearTimeout(onSlider.timer);
  onSlider.timer = setTimeout(() => {
    const image = state.images[state.selected];
    if (image) state.viewerImg.src = cropUrl(image, state.viewerRange);
  }, 120);
}

function refreshChannel(image) {
  for (const other of state.images) {
    if (other.pt === image.pt && other.channel === image.channel) other.img.src = cropUrl(other);
  }
  // the composite uses the channel's range too
  updateComposite();
}

// ---------- export ----------
function exportName(image) {
  // names the focused item: the cell and its well, or the dot and its tile and well
  const d = state.data;
  const tile = state.config.tiled ? `_${d.tile}` : '';
  const item = d.kind === 'cell' ? `${state.well}_cell${d.cell.id}` : `${state.well}${tile}_dot${d.index}`;
  const layers = compositeLayers();
  const what = image.composite
    ? (layers.length > 6 ? `composite${layers.length}` : 'composite_' + layers.map((l) => l.pt ? `${l.cycle}-${l.channel}` : `cycle${l.cycle}-${l.channel}`).join('+'))
    : image.pt ? `${image.cycle}_${image.channel}` : `cycle${image.cycle}_${image.channel}`;
  return `${item}_${what}`.replace(/[^\w.+-]+/g, '-');
}

const HATCH = ['#4a4a47', '#1f1f1e'];   // same as the page's 'not imaged' hatch

function overlaySvg(stage, size, withImage) {
  // one standalone SVG of the enlarged view: the annotation layers that are currently shown, with
  // colours and sizes resolved (CSS variables don't apply outside the page), scaled to size pixels
  const factor = size / stage.getBoundingClientRect().width;
  const root = svg('svg', { xmlns: SVG_NS, width: size, height: size, viewBox: `0 0 ${size} ${size}` });
  if (withImage) {
    const defs = svg('defs');
    const pattern = svg('pattern', { id: 'not-imaged', width: 6 * factor, height: 6 * factor, patternUnits: 'userSpaceOnUse', patternTransform: 'rotate(45)' });
    pattern.append(svg('rect', { width: 6 * factor, height: 6 * factor, fill: HATCH[1] }), svg('rect', { width: 2 * factor, height: 6 * factor, fill: HATCH[0] }));
    defs.append(pattern);
    root.append(defs, svg('rect', { width: size, height: size, fill: 'url(#not-imaged)' }));
    const image = svg('image', { href: withImage, width: size, height: size, preserveAspectRatio: 'none', 'image-rendering': 'pixelated' });
    root.append(image);
  }
  const shown = (node) => getComputedStyle(node).display !== 'none';
  for (const layer of stage.querySelectorAll(':scope > svg.overlay')) {
    if (!shown(layer)) continue;
    const copy = svg('svg', { x: 0, y: 0, width: size, height: size, viewBox: layer.getAttribute('viewBox'), preserveAspectRatio: 'none' });
    for (const node of layer.querySelectorAll('path, circle')) {
      if (node.classList.contains('region-hit')) continue;   // invisible click targets
      const style = getComputedStyle(node);
      const clone = node.cloneNode(false);
      clone.removeAttribute('class');
      clone.style.cssText = '';
      clone.setAttribute('stroke', style.stroke);
      clone.setAttribute('fill', 'none');
      clone.setAttribute('stroke-width', parseFloat(style.strokeWidth) * factor);
      copy.append(clone);
    }
    root.append(copy);
  }
  // dot labels and 'not imaged' text are HTML on the page; redraw them as SVG text in the same place
  const origin = stage.getBoundingClientRect();
  for (const text of stage.querySelectorAll('.dot-label, .no-data')) {
    if (!shown(text) || !text.textContent) continue;
    const style = getComputedStyle(text);
    const box = text.getBoundingClientRect();
    const fontSize = parseFloat(style.fontSize) * factor;
    const label = svg('text', {
      x: (box.left - origin.left) * factor + (text.classList.contains('no-data') ? box.width * factor / 2 : 0),
      y: (box.top - origin.top + box.height / 2) * factor,
      'font-family': style.fontFamily, 'font-size': fontSize, 'font-weight': style.fontWeight,
      'dominant-baseline': 'central', 'text-anchor': text.classList.contains('no-data') ? 'middle' : 'start',
      fill: style.color, stroke: '#000', 'stroke-width': fontSize / 5, 'paint-order': 'stroke', 'stroke-linejoin': 'round',
    });
    label.textContent = text.textContent;
    root.append(label);
  }
  return new XMLSerializer().serializeToString(root);
}

function cropDataUrl(img) {
  // the enlarged view's crop PNG at its own resolution
  const canvas = document.createElement('canvas');
  canvas.width = img.naturalWidth; canvas.height = img.naturalHeight;
  canvas.getContext('2d').drawImage(img, 0, 0);
  return canvas.toDataURL('image/png');
}

function download(blob, name) {
  const url = URL.createObjectURL(blob);
  const link = el('a', { href: url, download: name });
  document.body.append(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 10000);
}

async function exportView(format) {
  const image = state.images[state.selected];
  const stage = $('viewer-stage').firstElementChild;
  const img = state.viewerImg;
  if (!image || !stage || !img) return;
  const status = $('export-status');
  try {
    if (!img.complete || !img.naturalWidth) await new Promise((resolve, reject) => { img.onload = resolve; img.onerror = reject; });
    const name = exportName(image);
    if (format === 'svg') {
      const text = overlaySvg(stage, Math.round(stage.getBoundingClientRect().width), cropDataUrl(img));
      download(new Blob([text], { type: 'image/svg+xml' }), name + '.svg');
      status.textContent = `saved ${name}.svg`;
      return;
    }
    const size = parseInt($('export-size').value, 10);
    const canvas = document.createElement('canvas');
    canvas.width = canvas.height = size;
    const ctx = canvas.getContext('2d');
    // hatch under the crop so pixels that were never imaged read as missing, as on the page
    const tile = document.createElement('canvas');
    const step = Math.max(4, Math.round(6 * size / stage.getBoundingClientRect().width));
    tile.width = tile.height = step;
    const t = tile.getContext('2d');
    t.fillStyle = HATCH[1]; t.fillRect(0, 0, step, step);
    t.strokeStyle = HATCH[0]; t.lineWidth = step / 3;
    for (const offset of [-step, 0, step]) { t.beginPath(); t.moveTo(offset, step); t.lineTo(offset + step, 0); t.stroke(); }
    ctx.fillStyle = ctx.createPattern(tile, 'repeat'); ctx.fillRect(0, 0, size, size);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(img, 0, 0, size, size);
    const overlay = new Image();
    const url = URL.createObjectURL(new Blob([overlaySvg(stage, size, null)], { type: 'image/svg+xml' }));
    await new Promise((resolve, reject) => { overlay.onload = resolve; overlay.onerror = reject; overlay.src = url; });
    ctx.drawImage(overlay, 0, 0, size, size);
    URL.revokeObjectURL(url);
    canvas.toBlob((blob) => { download(blob, name + '.png'); status.textContent = `saved ${name}.png`; }, 'image/png');
  } catch (error) {
    status.textContent = 'export failed: ' + (error.message || error.type || error);
  }
}

// ---------- charts ----------
function niceTicks(lo, hi, count = 5) {
  const span = hi - lo || 1;
  const step0 = span / count;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => span / s <= count) || 10 * mag;
  const ticks = [];
  const last = Math.ceil(hi / step - 1e-9) * step;
  for (let v = Math.floor(lo / step + 1e-9) * step; v <= last + step * 1e-9; v += step) ticks.push(+v.toFixed(10));
  return ticks;
}

function lineChart({ title, cycles, series, domain, xLabels, digits = 3 }) {
  const W = 560, H = 300, m = { l: 50, r: series.length > 1 ? 30 : 14, t: 10, b: 46 };
  const box = el('div', { class: 'chart' }, el('h2', { text: title }));
  if (series.length > 1) {
    const legend = el('div', { class: 'legend' });
    for (const s of series) {
      const sw = el('i'); sw.style.background = s.color;
      legend.append(el('span', {}, sw, s.name));
    }
    box.append(legend);
  }
  const all = series.flatMap((s) => s.values).filter((v) => Number.isFinite(v));
  let [lo, hi] = domain || [Math.min(0, ...all), Math.max(0, ...all)];
  const ticks = niceTicks(lo, hi);
  lo = Math.min(lo, ticks[0]); hi = Math.max(hi, ticks[ticks.length - 1]);
  const x = (i) => m.l + (cycles.length === 1 ? 0.5 : i / (cycles.length - 1)) * (W - m.l - m.r);
  const y = (v) => m.t + (1 - (v - lo) / (hi - lo || 1)) * (H - m.t - m.b);

  const chart = svg('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': title });
  for (const t of ticks) {
    chart.append(svg('line', { class: t === 0 ? 'baseline' : 'gridline', x1: m.l, x2: W - m.r, y1: y(t), y2: y(t) }));
    const label = svg('text', { class: 'tick', x: m.l - 6, y: y(t) + 5, 'text-anchor': 'end' });
    label.textContent = t;
    chart.append(label);
  }
  cycles.forEach((c, i) => {
    const label = svg('text', { class: 'tick', x: x(i), y: H - m.b + 19, 'text-anchor': 'middle' });
    label.textContent = c;
    chart.append(label);
    if (xLabels && xLabels[i]) {
      const base = svg('text', { class: 'callbase', x: x(i), y: H - m.b + 38, 'text-anchor': 'middle' });
      base.textContent = xLabels[i];
      chart.append(base);
    }
  });
  const crosshair = svg('line', { class: 'crosshair', y1: m.t, y2: H - m.b, visibility: 'hidden' });
  chart.append(crosshair);

  const ends = [];
  for (const s of series) {
    const points = s.values.map((v, i) => [x(i), y(v)]).filter((p, i) => Number.isFinite(s.values[i]));
    chart.append(svg('path', { class: 'line', stroke: s.color, d: points.map((p, i) => (i ? 'L' : 'M') + p.join(',')).join('') }));
    for (const [px, py] of points) chart.append(svg('circle', { class: 'pt', cx: px, cy: py, r: 4, fill: s.color }));
    if (series.length > 1 && points.length) ends.push({ name: s.name, x: points[points.length - 1][0], y: points[points.length - 1][1] });
  }
  // direct labels at line ends, pushed apart so converging lines stay readable
  ends.sort((a, b) => a.y - b.y);
  for (let i = 1; i < ends.length; i++) ends[i].y = Math.max(ends[i].y, ends[i - 1].y + 15);
  for (const end of ends) {
    const label = svg('text', { class: 'direct', x: end.x + 8, y: end.y + 4 });
    label.textContent = end.name;
    chart.append(label);
  }

  // hover: one hit column per cycle, crosshair + tooltip with every series
  const tooltip = $('tooltip');
  const step = cycles.length > 1 ? (x(1) - x(0)) : (W - m.l - m.r);
  cycles.forEach((c, i) => {
    const hit = svg('rect', { class: 'hit', x: x(i) - step / 2, y: m.t, width: step, height: H - m.t - m.b });
    hit.addEventListener('mousemove', (e) => {
      crosshair.setAttribute('x1', x(i)); crosshair.setAttribute('x2', x(i));
      crosshair.setAttribute('visibility', 'visible');
      tooltip.replaceChildren(el('div', {}, el('b', { text: `cycle ${c}` }), xLabels && xLabels[i] ? ` · called ${xLabels[i]}` : ''));
      for (const s of series) {
        const sw = el('span', { class: 'sw' }); sw.style.background = s.color;
        tooltip.append(el('div', {}, sw, `${s.name}: ${fmt(s.values[i], digits)}`));
      }
      tooltip.hidden = false;
      tooltip.style.left = (e.clientX + 14) + 'px';
      tooltip.style.top = (e.clientY + 14) + 'px';
    });
    hit.addEventListener('mouseleave', () => { crosshair.setAttribute('visibility', 'hidden'); tooltip.hidden = true; });
    chart.append(hit);
  });
  box.append(chart);

  // table view of the same numbers
  const table = el('table', { class: 'data' });
  table.append(el('tr', {}, el('th', { text: 'cycle' }), ...series.map((s) => el('th', { text: s.name }))));
  cycles.forEach((c, i) => table.append(el('tr', {}, el('td', { text: c }), ...series.map((s) => el('td', { text: fmt(s.values[i], digits) })))));
  box.append(el('details', {}, el('summary', { text: 'Table' }), table));
  return box;
}

function renderDotCharts() {
  const d = state.data;
  const calls = [...d.sequence];
  $('charts').replaceChildren(
    lineChart({ title: 'Chastity per cycle', cycles: d.cycles, domain: [0, 1], xLabels: calls,
      series: [{ name: 'chastity', color: 'var(--series-1)', values: d.chastity }] }),
    lineChart({ title: 'Z-score per cycle', cycles: d.cycles, xLabels: calls, digits: 2,
      series: d.zscore_bases.map((b) => ({ name: b, color: baseColor(b), values: d.zscores[b] })) }),
  );
}

function renderCellExtra() {
  const d = state.data;
  const box = $('cell-extra');
  box.append(el('h2', { text: 'Reads' }));
  if (d.reads.length) {
    const table = el('table', { class: 'data' });
    table.append(el('tr', {}, ...['read', 'count', 'chastities', 'hamming', 'barcode matches'].map((t) => el('th', { text: t }))));
    for (const r of d.reads) {
      table.append(el('tr', {},
        el('td', { class: 'seq', text: r.read }), el('td', { text: r.count ? parseFloat(r.count) : '' }),
        el('td', { class: 'wrap', text: r.chastities || '' }), el('td', { text: r.barcode_hamming_dist ? parseFloat(r.barcode_hamming_dist) : '' }),
        el('td', { class: 'wrap seq', text: (r.barcode_matches || '').replaceAll(';', '; ') })));
    }
    box.append(table);
  } else {
    box.append(el('p', { class: 'note', text: 'No reads in this cell.' }));
  }

  const dotsHead = el('h2', { text: 'Dots in cell' });
  dotsHead.style.marginTop = '14px';   // a style property, not a style attribute, which the page's CSP blocks
  box.append(dotsHead);
  if (d.dots.length) {
    const table = el('table', { class: 'data' });
    table.append(el('tr', {}, ...['dot index', 'sequence', 'min chastity'].map((t) => el('th', { text: t }))));
    for (const dot of d.dots) {
      table.append(el('tr', { class: 'clickable', title: 'Open dot view', onclick: () => { location.hash = dotHash(d.tile, dot.index); } },
        el('td', { text: dot.index }), el('td', { class: 'seq', text: dot.sequence }), el('td', { text: fmt(dot.min_chastity) })));
    }
    box.append(table);
  } else {
    box.append(el('p', { class: 'note', text: 'No dots in this cell.' }));
  }
  box.append(el('p', { class: 'note', text: `Dots from ${d.dot_source}. Click a dot to open it.` }));
}

// ---------- wiring ----------
async function init() {
  try {
    state.config = await api('/api/config');
  } catch (error) {
    $('error').hidden = false;
    $('error').textContent = 'Could not reach the viewer server: ' + error.message;
    return;
  }
  const c = state.config;
  $('run-name').textContent = c.run_dir.split('/').filter(Boolean).pop();
  for (const well of c.wells) $('well').append(el('option', { value: well, text: well }));
  $('well').addEventListener('change', () => {
    // show the same dot or cell in the chosen well
    const target = parseHash(), well = $('well').value;
    if (target.mode === 'dot' && target.tile != null) location.hash = `#dot/${well}/${target.tile}/${target.index}` + (target.size ? `?size=${target.size}` : '');
    else if (target.mode === 'cell' && target.id != null) location.hash = `#cell/${well}/${target.id}` + (target.size ? `?size=${target.size}` : '');
    else location.hash = `#${state.mode}/${well}`;
  });

  $('tab-dot').addEventListener('click', () => setMode('dot'));
  $('tab-cell').addEventListener('click', () => setMode('cell'));
  $('dot-form').addEventListener('submit', (e) => {
    e.preventDefault();
    location.hash = dotHash($('dot-tile').value, $('dot-index').value, $('dot-size').value);
  });
  $('cell-form').addEventListener('submit', (e) => { e.preventDefault(); location.hash = cellHash($('cell-id').value, $('cell-size').value); });
  // row labels are measured, so measure again once Ubuntu Mono has loaded
  if (document.fonts) document.fonts.ready.then(fitRowHeads);
  $('contrast-mode').addEventListener('change', () => { if (state.data) refreshImages(); });
  // each option value lists what it shows ('cells-nuclei-dots-ids'); every part becomes a show-<part> class on <body>
  const overlayMode = $('overlay-mode');
  const OLD_MODES = { cell: 'cells', all: 'cells-dots', ids: 'cells-dots-ids', seqs: 'cells-dots-seqs' };
  let saved = null;
  try { saved = localStorage.getItem('starcall-viewer-overlays'); } catch (e) { /* storage unavailable */ }
  overlayMode.value = OLD_MODES[saved] || saved || 'cells-dots';
  if (!overlayMode.value) overlayMode.value = 'cells-dots';
  const applyOverlayMode = () => {
    document.body.classList.remove(...[...document.body.classList].filter((name) => name.startsWith('show-')));
    if (overlayMode.value !== 'none') document.body.classList.add(...overlayMode.value.split('-').map((part) => 'show-' + part));
    try { localStorage.setItem('starcall-viewer-overlays', overlayMode.value); } catch (e) { /* storage unavailable */ }
  };
  overlayMode.addEventListener('change', applyOverlayMode);
  $('show-region').checked = false;   // off on every page load
  $('show-region').addEventListener('change', (e) => {
    if (e.target.checked) loadRegion();
    else { state.region = null; $('region-status').textContent = ''; if (state.data) refreshOverlays(); }
  });
  applyOverlayMode();
  $('prev-image').addEventListener('click', () => select(state.selected - 1));
  $('next-image').addEventListener('click', () => select(state.selected + 1));
  $('vmin').addEventListener('input', onSlider);
  $('vmax').addEventListener('input', onSlider);
  $('apply-contrast').addEventListener('click', () => {
    const image = state.images[state.selected];
    if (!image) return;
    state.overrides[`${image.pt ? 'pt' : 'seq'}:${image.channel}`] = { vmin: parseFloat($('vmin').value), vmax: parseFloat($('vmax').value) };
    refreshChannel(image);
  });
  $('reset-colors').addEventListener('click', () => { state.colors = {}; saveColors(); applyColors(); });
  $('composite-default').addEventListener('click', () => {
    state.composite = null;
    saveComposite();
    if (state.data) { renderLayerPicker(); updateComposite(); }
  });
  $('export-png').addEventListener('click', () => exportView('png'));
  $('export-svg').addEventListener('click', () => exportView('svg'));
  $('reset-contrast').addEventListener('click', () => {
    const image = state.images[state.selected];
    if (!image) return;
    delete state.overrides[`${image.pt ? 'pt' : 'seq'}:${image.channel}`];
    refreshChannel(image);
    select(state.selected);
  });
  document.addEventListener('keydown', (e) => {
    if (e.target.matches('input, select')) return;
    if (e.key === 'ArrowLeft') select(state.selected - 1);
    if (e.key === 'ArrowRight') select(state.selected + 1);
  });
  window.addEventListener('hashchange', route);
  await route();
}

init();
