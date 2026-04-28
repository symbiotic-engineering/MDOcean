from pathlib import Path
from playwright.sync_api import sync_playwright
import subprocess, time

slides = Path(__file__).parent / "slides.html"
#url = slides.resolve().as_uri() + "?print-pdf&transition=none"
url = "http://localhost:8000/slides.html?print-pdf"
out = Path(__file__).parent / "slides.pdf"

server = subprocess.Popen(
    ["python", "-m", "http.server", "8000", "--directory", "/work/pubs/defense"]
)

time.sleep(2)

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    page.add_init_script("""
    window.RevealConfigure = { transition: 'none' };
    """)
    page.goto(url, wait_until="networkidle")

    # Wait for Reveal to exist
    page.wait_for_function("() => window.Reveal")

    # Wait for slide construction to finish
    page.wait_for_function("""
    () => {
    const r = window.Reveal;
    if (!r || !r.getSlides) return false;

    const slides = r.getSlides();
    return slides && slides.length > 0;
    }
    """)

    # Wait for DOM order stabilization
    page.wait_for_function("""
    () => {
    const sections = document.querySelectorAll('.reveal .slides section');
    return sections.length > 1;
    }
    """)

    slides = page.evaluate("""
    () => {
    const slides = Reveal.getSlides(); // IMPORTANT: Reveal canonical order

    return slides.map((s, i) => ({
        index: i,
        id: s.id || null,
        html: s.innerHTML,
        isTitle: s.id === "title-slide"
    }));
    }
    """)
    title = [s for s in slides if s["isTitle"]]
    rest = [s for s in slides if not s["isTitle"]]

    fixed = title + rest

    #fixed = slides[shift:shift+1] + slides[:shift] + slides[shift+1:]

    html = """
    <html>
    <head>
    <link rel="stylesheet" href="slides_files/libs/revealjs/dist/reveal.css">
    </head>
    <body>
    <div class="reveal">
    <div class="slides">
    """

    for s in fixed:
        html += f"<section>{s['html']}</section>"

    html += """
    </div>
    </div>
    </body>
    </html>
    """

    page.set_content(html, wait_until="load")
    page.wait_for_timeout(1000)

    # Final stabilization
    page.evaluate("""
    () => {
    Reveal.layout();
    Reveal.sync();
    Reveal.slide(0,0);
    }
    """)

    page.wait_for_timeout(1500)
    #page.emulate_media(media="screen")
    page.set_viewport_size({"width": 1920, "height": 1080})
    print(page.evaluate("""
    () => Array.from(document.querySelectorAll('.reveal .slides section'))
    .map(s => s.innerText.slice(0, 40))
    """))
    page.screenshot(path="debug.png", full_page=True)
    # small extra delay helps with fonts/layout
    page.pdf(
        path=str(out), 
        print_background=True,
        prefer_css_page_size=True
    )
    browser.close()

server.terminate()