// 3D union map: metabolites as points on three compartment shells, reactions as
// substrate to product edges, coloured by pathway group. three.js is loaded from
// the CDN importmap at runtime; when WebGL or the CDN is unavailable the caller
// renders a text fallback and the rest of the app keeps working.

export async function createMap(container, graph, groupColors, opts = {}) {
  const canvasTest = document.createElement('canvas');
  const gl = canvasTest.getContext('webgl2') || canvasTest.getContext('webgl');
  if (!gl) throw new Error('webgl-unavailable');

  let THREE, OrbitControls, CSS2DRenderer, CSS2DObject;
  try {
    THREE = await import('three');
    ({ OrbitControls } = await import('three/addons/controls/OrbitControls.js'));
    ({ CSS2DRenderer, CSS2DObject } = await import('three/addons/renderers/CSS2DRenderer.js'));
  } catch (e) {
    throw new Error('cdn-unavailable');
  }

  const reduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const W = () => container.clientWidth, H = () => container.clientHeight;

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0xffffff);
  const camera = new THREE.PerspectiveCamera(45, W() / H(), 0.1, 600);
  camera.position.set(62, 40, 78);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setSize(W(), H());
  container.appendChild(renderer.domElement);
  renderer.domElement.setAttribute('tabindex', '0');
  renderer.domElement.setAttribute('aria-label', '3D metabolic map. Drag to orbit, scroll to zoom.');

  const labelRenderer = new CSS2DRenderer();
  labelRenderer.setSize(W(), H());
  labelRenderer.domElement.style.position = 'absolute';
  labelRenderer.domElement.style.top = '0';
  labelRenderer.domElement.style.pointerEvents = 'none';
  container.appendChild(labelRenderer.domElement);

  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = !reduced;
  controls.dampingFactor = 0.08;
  controls.minDistance = 18;
  controls.maxDistance = 260;

  // ---- shells
  const shellDefs = [
    { r: graph.shells.e, label: `extracellular / exchange (r ${graph.shells.e})` },
    { r: graph.shells.p, label: `periplasm (r ${graph.shells.p})` },
    { r: graph.shells.c, label: `cytosol (r ${graph.shells.c})` },
  ];
  for (const s of shellDefs) {
    const geo = new THREE.SphereGeometry(s.r, 28, 18);
    const mat = new THREE.MeshBasicMaterial({ color: 0xE6E3DC, wireframe: true, transparent: true, opacity: 0.22 });
    scene.add(new THREE.Mesh(geo, mat));
    const div = document.createElement('div');
    div.className = 'shell-label';
    div.textContent = s.label;
    const lab = new CSS2DObject(div);
    lab.position.set(0, s.r + 2.5, 0);
    scene.add(lab);
  }

  // ---- node sprite texture (round points)
  const spriteCanvas = document.createElement('canvas');
  spriteCanvas.width = spriteCanvas.height = 64;
  const ctx = spriteCanvas.getContext('2d');
  ctx.beginPath(); ctx.arc(32, 32, 28, 0, Math.PI * 2);
  ctx.fillStyle = '#fff'; ctx.fill();
  const spriteTex = new THREE.CanvasTexture(spriteCanvas);

  // ---- nodes: main + currency
  const mids = Object.keys(graph.metabolites);
  const mainIds = [], curIds = [];
  for (const mid of mids) (graph.metabolites[mid].cur ? curIds : mainIds).push(mid);

  const col = new THREE.Color();
  function buildPoints(ids, size, opacity, dimGrey) {
    const pos = new Float32Array(ids.length * 3);
    const colors = new Float32Array(ids.length * 3);
    ids.forEach((mid, i) => {
      const m = graph.metabolites[mid];
      pos.set(m.p, i * 3);
      col.set(dimGrey ? '#B4B0A8' : (groupColors[m.g] || '#9AA0A6'));
      colors.set([col.r, col.g, col.b], i * 3);
    });
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    const mat = new THREE.PointsMaterial({
      size, vertexColors: true, map: spriteTex, alphaTest: 0.4,
      transparent: true, opacity, sizeAttenuation: true,
    });
    return new THREE.Points(geo, mat);
  }
  const mainPoints = buildPoints(mainIds, 1.1, 0.95, false);
  const curPoints = buildPoints(curIds, 0.55, 0.4, true);
  scene.add(mainPoints, curPoints);

  // ---- edges: substrate -> product segments; currency-touching pairs kept dim
  function buildEdges() {
    const mainPos = [], mainCol = [], curPos = [];
    for (const r of graph.reactions) {
      col.set(groupColors[r.g] || '#9AA0A6');
      for (const s of r.s) {
        const ms = graph.metabolites[s];
        if (!ms) continue;
        for (const p of r.p) {
          const mp = graph.metabolites[p];
          if (!mp) continue;
          if (ms.cur || mp.cur) {
            curPos.push(...ms.p, ...mp.p);
          } else {
            mainPos.push(...ms.p, ...mp.p);
            mainCol.push(col.r, col.g, col.b, col.r, col.g, col.b);
          }
        }
      }
    }
    const g1 = new THREE.BufferGeometry();
    g1.setAttribute('position', new THREE.BufferAttribute(new Float32Array(mainPos), 3));
    g1.setAttribute('color', new THREE.BufferAttribute(new Float32Array(mainCol), 3));
    const e1 = new THREE.LineSegments(g1, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.14 }));
    const g2 = new THREE.BufferGeometry();
    g2.setAttribute('position', new THREE.BufferAttribute(new Float32Array(curPos), 3));
    const e2 = new THREE.LineSegments(g2, new THREE.LineBasicMaterial({ color: 0xCCC8C0, transparent: true, opacity: 0.05 }));
    return [e1, e2];
  }
  const [mainEdges, curEdges] = buildEdges();
  scene.add(mainEdges, curEdges);

  // ---- pathway-group sector labels on the outer shell
  const sectorLabels = [];
  {
    const sums = new Map();
    for (const mid of mainIds) {
      const m = graph.metabolites[mid];
      const v = sums.get(m.g) || [0, 0, 0, 0];
      v[0] += m.p[0]; v[1] += m.p[1]; v[2] += m.p[2]; v[3]++;
      sums.set(m.g, v);
    }
    for (const g of graph.groups) {
      const v = sums.get(g);
      if (!v || v[3] < 3) continue;
      const len = Math.hypot(v[0], v[1], v[2]);
      if (len < 1e-6) continue;
      const R = graph.shells.e + 5;
      const div = document.createElement('div');
      div.className = 'sector-label';
      div.textContent = g;
      div.style.color = groupColors[g] ? '#1A1D21' : '';
      const lab = new CSS2DObject(div);
      lab.position.set(v[0] / len * R, v[1] / len * R, v[2] / len * R);
      scene.add(lab);
      sectorLabels.push(lab);
    }
  }

  // ---- hover tooltip + click focus
  const raycaster = new THREE.Raycaster();
  raycaster.params.Points = { threshold: 0.8 };
  const pointer = new THREE.Vector2();
  let hoverCb = opts.onHover || null;
  let pendingHover = null;

  function pick(ev) {
    const rect = renderer.domElement.getBoundingClientRect();
    pointer.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
    pointer.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(pointer, camera);
    const targets = curPoints.visible ? [mainPoints, curPoints] : [mainPoints];
    const hits = raycaster.intersectObjects(targets, false);
    if (!hits.length) return null;
    const h = hits[0];
    const mid = h.object === mainPoints ? mainIds[h.index] : curIds[h.index];
    return { mid, met: graph.metabolites[mid], x: ev.clientX - rect.left, y: ev.clientY - rect.top };
  }

  renderer.domElement.addEventListener('pointermove', (ev) => { pendingHover = ev; });
  let downAt = null;
  renderer.domElement.addEventListener('pointerdown', (ev) => { downAt = [ev.clientX, ev.clientY]; });
  renderer.domElement.addEventListener('pointerup', (ev) => {
    if (!downAt) return;
    const moved = Math.hypot(ev.clientX - downAt[0], ev.clientY - downAt[1]);
    downAt = null;
    if (moved > 4) return;
    const hit = pick(ev);
    if (hit) focusOn(hit.met.p);
  });

  let focusFrom = null, focusTo = null, focusT = 1;
  function focusOn(p) {
    focusFrom = controls.target.clone();
    focusTo = new THREE.Vector3(...p);
    focusT = reduced ? 1 : 0;
    if (reduced) controls.target.copy(focusTo);
  }

  // ---- pathway highlight overlay
  let hlGroup = null;
  function highlightPathway(metSeq, accentHex) {
    clearHighlight();
    hlGroup = new THREE.Group();
    const pts = metSeq.map(mid => graph.metabolites[mid]).filter(Boolean);
    if (pts.length < 2) return;
    const pos = [];
    for (let i = 0; i < pts.length - 1; i++) pos.push(...pts[i].p, ...pts[i + 1].p);
    const lg = new THREE.BufferGeometry();
    lg.setAttribute('position', new THREE.BufferAttribute(new Float32Array(pos), 3));
    const line = new THREE.LineSegments(lg, new THREE.LineBasicMaterial({ color: accentHex, transparent: true, opacity: 0.95, depthTest: false }));
    line.renderOrder = 10;
    const ng = new THREE.BufferGeometry();
    ng.setAttribute('position', new THREE.BufferAttribute(new Float32Array(pts.flatMap(m => m.p)), 3));
    const nodes = new THREE.Points(ng, new THREE.PointsMaterial({ color: accentHex, size: 2.2, map: spriteTex, alphaTest: 0.4, transparent: true, depthTest: false, sizeAttenuation: true }));
    nodes.renderOrder = 11;
    hlGroup.add(line, nodes);
    scene.add(hlGroup);
    mainEdges.material.opacity = 0.05;
    mainPoints.material.opacity = 0.35;
    curEdges.material.opacity = 0.02;
    focusOn(pts[Math.floor(pts.length / 2)].p);
  }
  function clearHighlight() {
    if (hlGroup) { scene.remove(hlGroup); hlGroup = null; }
    mainEdges.material.opacity = 0.14;
    mainPoints.material.opacity = 0.95;
    curEdges.material.opacity = 0.05;
  }

  // ---- render loop
  renderer.setAnimationLoop(() => {
    if (focusT < 1 && focusTo) {
      focusT = Math.min(1, focusT + 0.08);
      controls.target.lerpVectors(focusFrom, focusTo, 1 - Math.pow(1 - focusT, 3));
    }
    controls.update();
    if (pendingHover && hoverCb) {
      hoverCb(pick(pendingHover));
      pendingHover = null;
    }
    renderer.render(scene, camera);
    labelRenderer.render(scene, camera);
  });

  const ro = new ResizeObserver(() => {
    camera.aspect = W() / Math.max(H(), 1);
    camera.updateProjectionMatrix();
    renderer.setSize(W(), H());
    labelRenderer.setSize(W(), H());
  });
  ro.observe(container);

  return {
    highlightPathway,
    clearHighlight,
    setCurrencyVisible(v) { curPoints.visible = v; curEdges.visible = v; },
    setLabelsVisible(v) { sectorLabels.forEach(l => { l.visible = v; }); },
    resetView() {
      camera.position.set(62, 40, 78);
      controls.target.set(0, 0, 0);
      focusTo = null;
    },
    dispose() {
      renderer.setAnimationLoop(null);
      ro.disconnect();
      renderer.dispose();
      container.innerHTML = '';
    },
  };
}
