(() => {
  "use strict";
  const HPSC = window.HPSC || {};
  const modules = HPSC.modules || [];
  const index = HPSC.index || [];
  const slides = HPSC.slides || { decks: [] };
  const byId = Object.fromEntries(modules.map(m => [m.id, m]));
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => [...r.querySelectorAll(s)];
  const LS = "hpsc_progress_v1";
  const state = {
    current: null,
    progress: JSON.parse(localStorage.getItem(LS) || "{}"),
    deck: null,
    slideList: [],
    slidePos: 0,
  };

  /* ---------- helpers ---------- */
  function esc(s) { return String(s).replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c])); }
  function saveProgress() { localStorage.setItem(LS, JSON.stringify(state.progress)); }
  function moduleDone(id) {
    const m = byId[id]; if (!m) return false;
    const heads = m.headings || [];
    if (!heads.length) return !!state.progress[id];
    return heads.every(h => state.progress[id + "#" + h.id]) && !!state.progress[id];
  }
  function updateProgressUI() {
    const total = modules.length;
    const done = modules.filter(m => moduleDone(m.id)).length;
    const pct = total ? Math.round(100 * done / total) : 0;
    $("#progressbar").style.width = pct + "%";
    const ring = $("#doneCount"); if (ring) ring.textContent = done + "/" + total;
    $$(".navitem[data-module]").forEach(a => a.classList.toggle("is-done", moduleDone(a.dataset.module)));
    const cur = state.current;
    if (cur) {
      const heads = (byId[cur].headings || []);
      const cd = heads.filter(h => state.progress[cur + "#" + h.id]).length;
      const txt = $("#moduleProgressText");
      if (txt) txt.textContent = heads.length ? `${cd}/${heads.length} sections revised` : (state.progress[cur] ? "Revised" : "Not revised");
    }
  }

  /* ---------- navigation ---------- */
  function buildNav() {
    const nav = $("#nav");
    let html = "";
    let lastGroup = "";
    modules.forEach((m, i) => {
      if (m.group && m.group !== lastGroup) { html += `<div class="group">${esc(m.group)}</div>`; lastGroup = m.group; }
      html += `<a class="navitem" href="#${m.id}" data-module="${m.id}">
        <span class="num">${String(i).padStart(2, "0")}</span>
        <span>${esc(m.short || m.title)}</span><span class="done">✓</span></a>`;
    });
    nav.innerHTML = html;
    nav.addEventListener("click", e => {
      const a = e.target.closest("a[data-module]"); if (!a) return;
      if (window.innerWidth <= 860) document.body.classList.remove("nav-open");
    });
  }

  function updateSidebarSub(id) {
    $$(".sidebar .sub").forEach(s => s.classList.remove("open"));
    const m = byId[id]; if (!m) return;
    const sub = $("#sub-" + id); if (!sub) return;
    sub.classList.add("open");
    $$(".sidebar .sub a").forEach(a => a.classList.remove("active"));
    const first = sub.querySelector("a"); if (first) first.classList.add("active");
  }

  /* ---------- module rendering ---------- */
  function moduleSource(id) {
    const el = document.getElementById("src-" + id);
    return el ? el.textContent : "<p>Missing module.</p>";
  }

  function renderModule(id, opts = {}) {
    const m = byId[id];
    if (!m) { location.hash = modules[0].id; return; }
    const src = moduleSource(id);
    const content = $("#content");
    content.innerHTML = src;

    // wrap tables so wide ones scroll horizontally instead of shrinking
    $$("table", content).forEach(t => {
      const w = document.createElement("div"); w.className = "tablewrap";
      t.before(w); w.appendChild(t);
    });

    // header
    const head = document.createElement("header");
    head.className = "module-head";
    head.innerHTML = `<div class="eyebrow">${esc(m.group || "Module")} · ${esc(m.num)}</div>
      <h1>${esc(m.title)}</h1>
      ${m.lead ? `<p class="lead">${esc(m.lead)}</p>` : ""}
      <div class="module-meta">
        ${m.heat ? `<span class="badge ${m.heat === "Very high" ? "hot" : "ok"}">Exam probability: ${esc(m.heat)}</span>` : ""}
        <span class="badge">${m.mins || 10} min read</span>
        <span class="badge">${(m.words || 0).toLocaleString()} words</span>
        <label class="revised-toggle" style="margin-left:auto"><input type="checkbox" id="modDone" ${state.progress[id] ? "checked" : ""}> Module revised</label>
      </div>`;
    content.prepend(head);

    // toolbar + toc
    const heads = $$("h2", content).map((h, i) => {
      if (!h.id) h.id = id + "-h" + i;
      return { id: h.id, title: h.textContent.trim() };
    });
    m.headings = heads;
    const toc = document.createElement("nav");
    toc.className = "module-toc";
    toc.innerHTML = heads.map(h => `<a href="#${id}__${h.id}" data-sec="${h.id}">${esc(h.title)}</a>`).join("");
    const toolbar = document.createElement("div");
    toolbar.className = "toolbar";
    toolbar.innerHTML = `<button class="btn" id="printBtn">🖨 Print module</button>
      <button class="btn" id="topBtn">↑ Top</button>
      <span class="spacer"></span><span class="small" id="moduleProgressText"></span>`;
    head.after(toolbar); toolbar.after(toc);

    // per-section checkboxes
    heads.forEach((h, i) => {
      const el = document.getElementById(h.id); if (!el) return;
      const lab = document.createElement("label");
      lab.className = "revised-toggle";
      lab.style.cssText = "float:right;font-size:10px;margin-top:6px;opacity:.75";
      const cb = document.createElement("input"); cb.type = "checkbox"; cb.dataset.sec = h.id;
      cb.checked = !!state.progress[id + "#" + h.id];
      lab.appendChild(cb); lab.appendChild(document.createTextNode("revised"));
      el.appendChild(lab);
    });

    // copy buttons
    $$("pre", content).forEach(pre => {
      const b = document.createElement("button"); b.className = "copy"; b.textContent = "copy";
      b.onclick = () => { navigator.clipboard.writeText(pre.querySelector("code")?.innerText || pre.innerText).then(() => { b.textContent = "copied"; setTimeout(() => b.textContent = "copy", 1200); }); };
      pre.appendChild(b);
    });

    if (id === "slides") renderDeckList();

    // lazy images
    $$("img[data-src]", content).forEach(img => {
      if ("IntersectionObserver" in window) {
        const io = new IntersectionObserver((ents, obs) => ents.forEach(en => { if (en.isIntersecting) { img.src = img.dataset.src; img.removeAttribute("data-src"); obs.unobserve(img); } }), { rootMargin: "500px" });
        io.observe(img);
      } else img.src = img.dataset.src;
    });

    // wire toolbar
    $("#printBtn").onclick = () => window.print();
    $("#topBtn").onclick = () => window.scrollTo({ top: 0, behavior: "smooth" });
    $("#modDone").onchange = e => {
      state.progress[id] = e.target.checked;
      heads.forEach(h => { if (e.target.checked) state.progress[id + "#" + h.id] = true; else delete state.progress[id + "#" + h.id]; });
      $$("input[data-sec]", content).forEach(cb => cb.checked = e.target.checked);
      saveProgress(); updateProgressUI();
    };
    // section checkboxes are handled by the delegated listener in init()

    // toc + section tracking
    toc.addEventListener("click", e => {
      const a = e.target.closest("a[data-sec]"); if (!a) return;
      e.preventDefault();
      const el = document.getElementById(a.dataset.sec);
      if (el) window.scrollTo({ top: el.getBoundingClientRect().top + window.scrollY - 72, behavior: "smooth" });
    });

    state.current = id;
    $$(".navitem[data-module]").forEach(a => a.classList.toggle("active", a.dataset.module === id));
    updateSidebarSub(id);
    updateProgressUI();
    document.title = m.title + " · HPSC End-Sem Study Hub";
    if (opts.scroll !== false) window.scrollTo({ top: 0, behavior: "instant" in window ? "instant" : "auto" });
    if (!opts.keepHash && location.hash !== "#" + id) history.replaceState(null, "", "#" + id);
    localStorage.setItem("hpsc_last", id);
  }

  function route() {
    const full = (location.hash || "").replace(/^#/, "");
    const [id, sec] = full.split("__");
    const target = byId[id] ? id : modules[0].id;
    if (target !== state.current) renderModule(target, { scroll: !sec, keepHash: !!sec });
    function scrollToSection() {
      if (!sec) return;
      const el = document.getElementById(sec);
      if (!el) return;
      const y = el.getBoundingClientRect().top + window.scrollY - 72;
      const html = document.documentElement;
      const prev = html.style.scrollBehavior;
      html.style.scrollBehavior = "auto";
      window.scrollTo(0, y);
      html.scrollTop = y;
      document.body.scrollTop = y;
      html.style.scrollBehavior = prev;
    }
    // Scroll now, then again after MathML/images settle. Without the second
    // pass, a deep link can land above the heading while the page is still
    // reflowing (especially on a cold load).
    scrollToSection();
    if (sec) {
      requestAnimationFrame(scrollToSection);
      setTimeout(scrollToSection, 180);
    }
    const el = sec ? document.getElementById(sec) : null;
    if (location.search.includes("debug")) {
      const d = document.createElement("div");
      d.style.cssText = "position:fixed;left:10px;bottom:10px;z-index:999;background:#000;color:#0f0;font:12px monospace;padding:8px;max-width:90vw";
      d.textContent = `hash=${location.hash} | id=${id} | sec=${sec} | el=${el ? el.offsetTop : "NULL"} | scrollY=${window.scrollY}`;
      document.body.appendChild(d);
    }
  }

  /* ---------- search ---------- */
  let searchSel = 0, searchHits = [];
  function openSearch() { $("#searchmodal").classList.add("open"); const i = $("#searchInput"); i.value = ""; i.focus(); runSearch(""); }
  function closeSearch() { $("#searchmodal").classList.remove("open"); }
  function runSearch(q) {
    q = q.trim().toLowerCase();
    const toks = q.split(/\s+/).filter(Boolean);
    let hits = index;
    if (toks.length) hits = index.map(e => {
      const hay = (e.t + " " + e.s + " " + (e.k || "")).toLowerCase();
      let score = 0; toks.forEach(t => { if (e.t.toLowerCase().includes(t)) score += 6; else if (hay.includes(t)) score += 2; });
      return { e, score };
    }).filter(x => x.score > 0).sort((a, b) => b.score - a.score).map(x => x.e).slice(0, 40);
    else hits = index.slice(0, 24);
    searchHits = hits; searchSel = 0;
    const box = $("#searchResults");
    if (!hits.length) { box.innerHTML = `<div class="noresults">No match. Try "warp divergence", "halo", "Amdahl", "atomic".</div>`; return; }
    box.innerHTML = hits.map((e, i) => `<div class="sr ${i === 0 ? "sel" : ""}" data-i="${i}">
      <b>${esc(e.t)}</b><span>${esc(byId[e.m]?.short || e.m)}</span><p>${esc(e.s)}</p></div>`).join("");
  }
  function gotoHit(i) {
    const e = searchHits[i]; if (!e) return;
    closeSearch();
    const target = "#" + e.m + (e.id ? "__" + e.id : "");
    if (state.current === e.m) {
      const el = e.id ? document.getElementById(e.id) : null;
      if (el) window.scrollTo({ top: el.getBoundingClientRect().top + window.scrollY - 72, behavior: "smooth" });
      history.replaceState(null, "", target);
    } else location.hash = target;
  }

  /* ---------- slides ---------- */
  function renderDeckList() {
    const wrap = document.getElementById("deckList"); if (!wrap) return;
    const groups = {};
    slides.decks.forEach(d => { (groups[d.group] = groups[d.group] || []).push(d); });
    wrap.innerHTML = Object.entries(groups).map(([g, ds]) => `<div class="group" style="margin:14px 0 6px;font-size:10.5px;letter-spacing:1.2px;text-transform:uppercase;color:var(--muted);font-weight:800">${esc(g)}</div>
      <div class="deckgrid">${ds.map(d => `<div class="deck" data-deck="${esc(d.id)}"><b>${esc(d.name)}</b><span>${d.count} pages</span></div>`).join("")}</div>`).join("");
    wrap.addEventListener("click", e => {
      const el = e.target.closest(".deck"); if (!el) return;
      $$(".deck").forEach(x => x.classList.remove("active")); el.classList.add("active");
      showDeck(el.dataset.deck);
    });
  }
  function showDeck(id) {
    const d = slides.decks.find(x => x.id === id); if (!d) return;
    state.deck = id; state.slideList = d.pages || []; state.slidePos = 0;
    const grid = document.getElementById("slideGrid");
    grid.innerHTML = state.slideList.map((p, i) => `<div class="slidecard" data-i="${i}"><img src="assets/slides/${id}/${p}" loading="lazy" alt="${esc(d.name)} page ${i + 1}" width="320" height="200"><span class="cap">p.${i + 1}</span></div>`).join("");
    document.getElementById("deckTitle").textContent = d.name;
    document.getElementById("deckCount").textContent = ` · ${d.count} pages`;
    grid.onclick = e => { const c = e.target.closest(".slidecard"); if (c) openLightbox(+c.dataset.i); };
  }
  function openLightbox(i) {
    if (!state.slideList.length) return;
    state.slidePos = (i + state.slideList.length) % state.slideList.length;
    const d = slides.decks.find(x => x.id === state.deck);
    const img = document.getElementById("lbImg");
    img.src = `assets/slides/${state.deck}/${state.slideList[state.slidePos]}`;
    img.alt = d.name + " page " + (state.slidePos + 1);
    document.getElementById("lbCap").textContent = `${d.name} — page ${state.slidePos + 1} / ${d.count}`;
    document.getElementById("lbOpen").href = img.src;
    document.getElementById("lightbox").classList.add("open");
  }
  function closeLightbox() { document.getElementById("lightbox").classList.remove("open"); }

  /* ---------- theme ---------- */
  function setTheme(t) { document.documentElement.dataset.theme = t; localStorage.setItem("hpsc_theme", t); $("#themeBtn").textContent = t === "dark" ? "☀" : "☾"; }
  function toggleTheme() { setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark"); }

  /* ---------- init ---------- */
  function init() {
    buildNav();
    renderDeckList();
    setTheme(localStorage.getItem("hpsc_theme") || "dark");

    $("#searchBtn").onclick = openSearch;
    $("#searchInput").addEventListener("input", e => runSearch(e.target.value));
    $("#searchInput").addEventListener("keydown", e => {
      if (e.key === "ArrowDown") { e.preventDefault(); searchSel = Math.min(searchSel + 1, searchHits.length - 1); paintSel(); }
      if (e.key === "ArrowUp") { e.preventDefault(); searchSel = Math.max(searchSel - 1, 0); paintSel(); }
      if (e.key === "Enter") gotoHit(searchSel);
    });
    $("#searchResults").addEventListener("click", e => { const r = e.target.closest(".sr"); if (r) gotoHit(+r.dataset.i); });
    $("#searchmodal").addEventListener("click", e => { if (e.target.id === "searchmodal") closeSearch(); });
    function paintSel() { $$(".sr").forEach((r, i) => r.classList.toggle("sel", i === searchSel)); const s = $(".sr.sel"); if (s) s.scrollIntoView({ block: "nearest" }); }

    $("#themeBtn").onclick = toggleTheme;
    $("#content").addEventListener("change", e => {
      const cb = e.target.closest("input[data-sec]");
      if (!cb || !state.current) return;
      const key = state.current + "#" + cb.dataset.sec;
      if (cb.checked) state.progress[key] = true; else delete state.progress[key];
      saveProgress(); updateProgressUI();
    });
    $("#menuToggle").onclick = () => document.body.classList.toggle("nav-open");
    document.addEventListener("click", e => { if (e.target.classList.contains("mobile-scrim")) document.body.classList.remove("nav-open"); });
    window.addEventListener("keydown", e => {
      if (e.key === "/" && !/input|textarea/i.test(document.activeElement.tagName)) { e.preventDefault(); openSearch(); }
      if (e.key === "Escape") { closeSearch(); closeLightbox(); document.body.classList.remove("nav-open"); }
      if ($("#lightbox").classList.contains("open")) {
        if (e.key === "ArrowRight") openLightbox(state.slidePos + 1);
        if (e.key === "ArrowLeft") openLightbox(state.slidePos - 1);
      }
    });
    $("#lightbox").addEventListener("click", e => { if (e.target.id === "lightbox") closeLightbox(); });
    $("#lbPrev").onclick = () => openLightbox(state.slidePos - 1);
    $("#lbNext").onclick = () => openLightbox(state.slidePos + 1);
    $("#lbClose").onclick = closeLightbox;

    window.addEventListener("hashchange", route);
    window.addEventListener("load", () => { if (location.hash.includes("__")) route(); });
    window.addEventListener("scroll", () => {
      const h = document.documentElement;
      const p = h.scrollTop / Math.max(1, h.scrollHeight - h.clientHeight);
      if (!state.current) $("#progressbar").style.width = (p * 100).toFixed(1) + "%";
    }, { passive: true });

    const last = localStorage.getItem("hpsc_last");
    if (!location.hash && last && byId[last]) location.hash = last;
    route();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init); else init();
})();
