from pathlib import Path
from playwright.sync_api import sync_playwright
import subprocess, time

out = Path(__file__).parent / "slides.pdf"
slides = Path(__file__).parent / "slides.html"
url = slides.resolve().as_uri() + "?print-pdf"

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()

    page.goto(url, wait_until="networkidle")
    page.set_viewport_size({"width": 1280, "height": 720})
    page.wait_for_function("() => window.Reveal ")
    page.wait_for_function("() => document.querySelector('.reveal')")
    page.evaluate("""
    () => {
        document.body.classList.add('print-pdf');
        Reveal.layout();
        Reveal.sync();
    }
    """)

    # Force Reveal to recompute layout (THIS FIXES ORDER BUG)
    # page.evaluate("""
    # () => new Promise(resolve => {
    #     Reveal.on('ready', () => {
    #         Reveal.layout();
    #         Reveal.sync();
    #         resolve();
    #     });
    # })
    # """)

    page.emulate_media(media="print")

    # Give time for final reflow (important)
    page.wait_for_timeout(1500)
    page.pause()

    # debug sanity check
    headers = page.evaluate("""
    () => {
        return Array.from(document.querySelectorAll('.reveal h1, .reveal h2, .reveal h3'))
            .slice(0, 5)
            .map(el => {
                const style = window.getComputedStyle(el);
                return {
                    text: el.innerText,
                    color: style.color,
                    fontSize: style.fontSize,
                    parent: el.parentElement?.className
                };
            });
    }
    """)

    print(headers)

    debug = page.evaluate("""
    () => {
        const el = document.querySelector('.reveal h1');
        if (!el) return null;

        const style = getComputedStyle(el);

        return {
            color: style.color,
            cssText: style.cssText,
            className: el.className,
            parentClass: el.parentElement?.className,
            dataset: el.dataset
        };
    }
    """)

    print(debug)

    page.screenshot(path="debug.png", full_page=True)
    page.pdf(
        path=str(out),
        print_background=True,
        prefer_css_page_size=False,
        width="10in",
        height="5.625in",
        page_ranges="4-999"
    )

    browser.close()