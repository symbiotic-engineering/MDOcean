from pathlib import Path
from playwright.sync_api import sync_playwright

out = Path(__file__).parent / "slides.pdf"
slides = Path(__file__).parent / "slides.html"
url = slides.resolve().as_uri() + "?print-pdf"

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()

    page.goto(url, wait_until="networkidle")
    page.wait_for_function("() => window.Reveal ")
    page.wait_for_function("() => document.querySelector('.reveal')")
    page.evaluate("""
    () => {
        document.body.classList.add('print-pdf');
        Reveal.layout();
        Reveal.sync();
    }
    """)

    page.emulate_media(media="print")

    # Give time for final reflow (important)
    page.wait_for_timeout(1500)

    page.pdf(
        path=str(out),
        print_background=True,
        page_ranges="4-999"
    )

    browser.close()