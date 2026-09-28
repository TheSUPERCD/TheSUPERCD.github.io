#!/usr/bin/env python3
"""
update_site_nav.py
Updates INDEX.html, chapter-19.html, and likely-exam-questions.html to integrate Chapter 20.
"""
import os

def update_index_html(base_dir):
    index_file = os.path.join(base_dir, "INDEX.html")
    with open(index_file, "r", encoding="utf-8") as f:
        content = f.read()

    # 1. Update quick jump select
    old_opt = '<option value="chapter-19">Chapter 19: The Ultimate Last-Night Cramming Guide &amp; Complete Notation Decryptor</option>'
    new_opt = (
        '<option value="chapter-19">Chapter 19: The Ultimate Last-Night Cramming Guide &amp; Complete Notation Decryptor</option>\n'
        '<option value="chapter-20">Chapter 20: Sparse Matrix Storage Formats &amp; High-Performance Kernels (CSR, CSC, DIA, ELLPACK, ELLPACK-ITPACK, COO)</option>'
    )
    if old_opt in content and 'value="chapter-20"' not in content:
        content = content.replace(old_opt, new_opt, 1)
        print("Updated quick-jump select in INDEX.html")

    # 2. Update portal TOC list
    old_portal_li = '<li><a href="#chapter-19" onclick="closeTocMobile()">⚡ Chapter 19: Last-Night Cramming Guide &amp; Notation Decryptor</a></li>'
    new_portal_li = (
        '<li><a href="#chapter-19" onclick="closeTocMobile()">⚡ Chapter 19: Last-Night Cramming Guide &amp; Notation Decryptor</a></li>\n'
        '<li><a href="#chapter-20" onclick="closeTocMobile()">💾 Chapter 20: Sparse Storage Schemes &amp; Hardware Kernels</a></li>'
    )
    if old_portal_li in content and 'href="#chapter-20"' not in content:
        content = content.replace(old_portal_li, new_portal_li, 1)
        print("Updated portal TOC links in INDEX.html")

    # 3. Update TOC Sidebar: Add chapter-20 toc-group
    toc_marker = '<div class="toc-group" data-chapter="likely-exam-questions"'
    chapter_20_toc = """<div class="toc-group" data-chapter="chapter-20" data-chapter-id="chapter-20" data-search-text="Chapter 20 Sparse Matrix Storage Formats High-Performance Kernels CSR CSC DIA ELLPACK ELLPACK-ITPACK Ellpack-ltpack COO MSR 1. Foundations of Sparse Representations in Scientific Computing 2. Coordinate Format (COO / Triplet Format) 3. Compressed Sparse Row (CSR / CRS / Yale Format) 4. Compressed Sparse Column (CSC / CCS / Harwell-Boeing Format) 5. Diagonal Storage Format (DIA / DIAG) 6. ELLPACK (ELL Format) 7. ELLPACK-ITPACK Format (and ELLPACK-LTPACK Disambiguation) 8. Modified Sparse Row (MSR) Format 9. Comprehensive Comparative Analysis & Memory Engineering 10. Production-Grade Conversion Algorithms 11. High-Yield Examination Problems & Model Solutions">
<div class="toc-group-header">
<button aria-label="Toggle subtopics" class="toc-accordion-btn" onclick="toggleTocSub('toc-sub-chapter-20', this)" type="button">▸</button>
<a class="toc-group-title" href="#chapter-20" onclick="closeTocMobile()">
<span class="toc-badge badge-primary">Ch 20</span>
<span class="toc-title-text">Sparse Storage Schemes: CSR, CSC, DIA, ELLPACK, ELLPACK-ITPACK, COO</span>
</a>
<a class="toc-standalone-link" href="chapters/chapter-20.html" rel="noopener" target="_blank" title="Open standalone page for chapter-20">↗</a></div>
<ul class="toc-sub-list" id="toc-sub-chapter-20">
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#1-foundations-of-sparse-representations-in-scientific-computing" onclick="closeTocMobile()">1. Foundations of Sparse Representations</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#2-coordinate-format-coo--triplet-format" onclick="closeTocMobile()">2. Coordinate Format (COO / Triplet Format)</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#3-compressed-sparse-row-csr--crs--yale-format" onclick="closeTocMobile()">3. Compressed Sparse Row (CSR / CRS)</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#4-compressed-sparse-column-csc--ccs--harwell-boeing-format" onclick="closeTocMobile()">4. Compressed Sparse Column (CSC / CCS)</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#5-diagonal-storage-format-dia--diag" onclick="closeTocMobile()">5. Diagonal Storage Format (DIA / DIAG)</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#6-ellpack-ell-format" onclick="closeTocMobile()">6. ELLPACK (ELL Format)</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#7-ellpack-itpack-format-and-ellpack-ltpack-disambiguation" onclick="closeTocMobile()">7. ELLPACK-ITPACK Format (Ellpack-ltpack)</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#8-modified-sparse-row-msr-format" onclick="closeTocMobile()">8. Modified Sparse Row (MSR) Format</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#9-comprehensive-comparative-analysis--memory-engineering" onclick="closeTocMobile()">9. Comparative Analysis &amp; Memory Engineering</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#10-production-grade-conversion-algorithms" onclick="closeTocMobile()">10. Production-Grade Conversion Algorithms</a></li>
<li class="toc-depth-2"><a class="toc-target-link" data-chapter="chapter-20" href="#11-high-yield-examination-problems--model-solutions" onclick="closeTocMobile()">11. High-Yield Examination Problems</a></li>
</ul>
</div>
"""
    if toc_marker in content and 'data-chapter="chapter-20"' not in content:
        content = content.replace(toc_marker, chapter_20_toc + toc_marker, 1)
        print("Updated TOC Sidebar in INDEX.html")

    # 4. Update Dashboard Overview Table
    table_marker = "</tbody>\n</table>"
    ch20_row = """<tr>
<td style="text-align: left;"><a href="#chapter-20"><strong>Chapter 20</strong></a></td>
<td style="text-align: left;"><strong>Sparse Matrix Storage Formats &amp; High-Performance Computing Kernels</strong><br/>• Mathematical foundations of matrix sparsity; sparsity ratio (<span class="math inline">\(S\)</span>) vs. density (<span class="math inline">\(\\rho\)</span>).<br/>• The 6 Canonical Formats: Coordinate Format (COO), Compressed Sparse Row (CSR), Compressed Sparse Column (CSC), Diagonal Storage (DIA), ELLPACK, ELLPACK-ITPACK (Ellpack-ltpack), and Modified Sparse Row (MSR).<br/>• Detailed worked examples on canonical unstructured and structured matrices.<br/>• High-performance SpMV algorithms in C99, OpenMP parallelization, and GPU memory coalescing.<br/>• Rigorous memory footprint formulas, Roofline arithmetic intensity (<span class="math inline">\(I \\approx 0.16\)</span> FLOP/B), and breakeven density analysis.<br/>• Master Sparse Matrix Rosetta Stone &amp; Architecture-Aware Format Selection Decision Tree.<br/>• Complete model examination problems and solutions.</td>
<td style="text-align: center;"><span class="hub-badge badge-primary">Deep Dive</span><br/>Module 1 &amp; Hardware Kernels</td>
</tr>
"""
    if 'href="#chapter-20"><strong>Chapter 20</strong>' not in content and table_marker in content:
        content = content.replace(table_marker, ch20_row + table_marker, 1)
        print("Updated Dashboard Table in INDEX.html")

    # 5. Update Chapter Cards container
    card_marker = '<article class="chapter-card" data-chapter="Question Bank" id="likely-exam-questions">'
    ch20_card = """<article class="chapter-card" data-chapter="Chapter 20" id="chapter-20">
<header class="chapter-card-header">
<div class="chapter-header-left">
<div class="chapter-badge-row">
<span class="hub-badge badge-primary">Chapter 20</span>
<span class="hub-badge badge-part">Part I / Hardware Kernels</span>
</div>
<h2 class="chapter-title-heading">Sparse Matrix Storage Formats &amp; High-Performance Kernels: CSR, CSC, DIA, ELLPACK, ELLPACK-ITPACK, and COO</h2>
</div>
<div class="chapter-actions">
<button class="btn-action chapter-toggle-btn" onclick="toggleChapter('chapter-20')" title="Collapse or Expand this chapter" type="button"><span class="toggle-icon">＋</span> <span class="toggle-text">Expand</span></button><a class="btn-action btn-standalone" href="chapters/chapter-20.html" rel="noopener" target="_blank" title="Open in dedicated standalone page"><span>↗ Standalone</span></a>
<a class="btn-action btn-toc" href="#toc-sidebar" onclick="openToc()" title="Jump to Table of Contents">📑 TOC</a>
<a class="btn-action btn-top" href="#top" title="Scroll to top of document">⬆ Top</a>
</div>
</header>
<div class="chapter-card-body collapsed" data-loaded="false" id="chapter-20-content"></div>
</article>
"""
    if card_marker in content and 'id="chapter-20"' not in content:
        content = content.replace(card_marker, ch20_card + card_marker, 1)
        print("Updated Chapter Cards in INDEX.html")

    with open(index_file, "w", encoding="utf-8") as f:
        f.write(content)
    print("Saved INDEX.html")


def update_chapter_19(base_dir):
    ch19_file = os.path.join(base_dir, "chapters", "chapter-19.html")
    with open(ch19_file, "r", encoding="utf-8") as f:
        content = f.read()

    # Top nav Next link
    old_top_next = '<a href="likely-exam-questions.html" class="btn-nav-control" title="Next: HPSC Examination Question Bank & Model Solutions Repository">Next →</a>'
    new_top_next = '<a href="chapter-20.html" class="btn-nav-control" title="Next: Chapter 20: Sparse Matrix Storage Formats & High-Performance Kernels">Next →</a>'
    if old_top_next in content:
        content = content.replace(old_top_next, new_top_next, 1)
        print("Updated top navbar in chapter-19.html")

    # Bottom nav card
    old_bottom_next = """      <a href="likely-exam-questions.html" class="nav-card nav-next">
        <span class="nav-card-label">Next Resource →</span>
        <span class="nav-card-title">HPSC Examination Question Bank & Model Solutions Repository</span>
      </a>"""
    new_bottom_next = """      <a href="chapter-20.html" class="nav-card nav-next">
        <span class="nav-card-label">Next Chapter →</span>
        <span class="nav-card-title">Chapter 20: Sparse Matrix Storage Formats & High-Performance Kernels (CSR, CSC, DIA, ELLPACK, ELLPACK-ITPACK, COO)</span>
      </a>"""
    if old_bottom_next in content:
        content = content.replace(old_bottom_next, new_bottom_next, 1)
        print("Updated bottom nav card in chapter-19.html")

    with open(ch19_file, "w", encoding="utf-8") as f:
        f.write(content)
    print("Saved chapter-19.html")


def update_likely_exam_questions(base_dir):
    leq_file = os.path.join(base_dir, "chapters", "likely-exam-questions.html")
    with open(leq_file, "r", encoding="utf-8") as f:
        content = f.read()

    # Top nav Prev link
    old_top_prev = '<a href="chapter-19.html" class="btn-nav-control" title="Previous: The Ultimate Last-Night Cramming Guide & Complete Notation Decryptor">← Prev</a>'
    new_top_prev = '<a href="chapter-20.html" class="btn-nav-control" title="Previous: Chapter 20: Sparse Matrix Storage Formats & High-Performance Kernels">← Prev</a>'
    if old_top_prev in content:
        content = content.replace(old_top_prev, new_top_prev, 1)
        print("Updated top navbar in likely-exam-questions.html")

    # Bottom nav card
    old_bottom_prev = """      <a href="chapter-19.html" class="nav-card nav-prev">
        <span class="nav-card-label">← Previous Chapter</span>
        <span class="nav-card-title">Chapter 19: The Ultimate Last-Night Cramming Guide & Complete Notation Decryptor</span>
      </a>"""
    new_bottom_prev = """      <a href="chapter-20.html" class="nav-card nav-prev">
        <span class="nav-card-label">← Previous Chapter</span>
        <span class="nav-card-title">Chapter 20: Sparse Matrix Storage Formats & High-Performance Kernels</span>
      </a>"""
    if old_bottom_prev in content:
        content = content.replace(old_bottom_prev, new_bottom_prev, 1)
        print("Updated bottom nav card in likely-exam-questions.html")

    with open(leq_file, "w", encoding="utf-8") as f:
        f.write(content)
    print("Saved likely-exam-questions.html")


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    update_index_html(base_dir)
    update_chapter_19(base_dir)
    update_likely_exam_questions(base_dir)

if __name__ == "__main__":
    main()
