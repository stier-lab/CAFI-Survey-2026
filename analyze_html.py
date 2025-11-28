#!/usr/bin/env python3
"""
Analyze and critique HTML output from RMarkdown files using Playwright
"""

import asyncio
from pathlib import Path
from playwright.async_api import async_playwright
import json

async def analyze_html(html_path):
    """Analyze an HTML file for design and usability issues"""

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()

        # Load the HTML file
        file_url = f"file://{Path(html_path).absolute()}"
        await page.goto(file_url)

        # Wait for content to load
        await page.wait_for_load_state('networkidle')

        # Analyze various aspects
        analysis = {}

        # 1. Check viewport and responsiveness
        viewport = page.viewport_size
        analysis['viewport'] = viewport

        # 2. Analyze typography
        analysis['typography'] = await page.evaluate('''() => {
            const body = document.body;
            const computedStyle = window.getComputedStyle(body);
            const headings = document.querySelectorAll('h1, h2, h3, h4, h5, h6');
            const paragraphs = document.querySelectorAll('p');

            return {
                bodyFontSize: computedStyle.fontSize,
                bodyLineHeight: computedStyle.lineHeight,
                bodyFontFamily: computedStyle.fontFamily,
                headingCount: headings.length,
                paragraphCount: paragraphs.length,
                avgParagraphLength: Array.from(paragraphs).reduce((acc, p) => acc + p.textContent.length, 0) / paragraphs.length
            };
        }''')

        # 3. Analyze images and figures
        analysis['images'] = await page.evaluate('''() => {
            const images = document.querySelectorAll('img');
            return {
                count: images.length,
                sizes: Array.from(images).map(img => ({
                    naturalWidth: img.naturalWidth,
                    naturalHeight: img.naturalHeight,
                    displayWidth: img.clientWidth,
                    displayHeight: img.clientHeight,
                    hasAlt: !!img.alt,
                    src: img.src.substring(img.src.lastIndexOf('/') + 1)
                }))
            };
        }''')

        # 4. Analyze navigation and structure
        analysis['navigation'] = await page.evaluate('''() => {
            const toc = document.querySelector('.toc, #TOC, [class*="toc"]');
            const sections = document.querySelectorAll('section, .section');
            const anchors = document.querySelectorAll('a[href^="#"]');

            return {
                hasTOC: !!toc,
                sectionCount: sections.length,
                internalLinkCount: anchors.length
            };
        }''')

        # 5. Analyze color scheme
        analysis['colors'] = await page.evaluate('''() => {
            const elements = document.querySelectorAll('*');
            const colors = new Set();
            const bgColors = new Set();

            elements.forEach(el => {
                const style = window.getComputedStyle(el);
                if (style.color) colors.add(style.color);
                if (style.backgroundColor && style.backgroundColor !== 'rgba(0, 0, 0, 0)') {
                    bgColors.add(style.backgroundColor);
                }
            });

            return {
                textColors: Array.from(colors).slice(0, 10),
                backgroundColors: Array.from(bgColors).slice(0, 10)
            };
        }''')

        # 6. Check for interactive elements
        analysis['interactivity'] = await page.evaluate('''() => {
            const buttons = document.querySelectorAll('button');
            const inputs = document.querySelectorAll('input, select, textarea');
            const details = document.querySelectorAll('details');
            const codeBlocks = document.querySelectorAll('pre, code');

            return {
                buttonCount: buttons.length,
                formElementCount: inputs.length,
                collapsibleCount: details.length,
                codeBlockCount: codeBlocks.length
            };
        }''')

        # 7. Measure loading performance
        analysis['performance'] = await page.evaluate('''() => {
            const perf = performance.timing;
            return {
                domContentLoaded: perf.domContentLoadedEventEnd - perf.navigationStart,
                loadComplete: perf.loadEventEnd - perf.navigationStart
            };
        }''')

        # 8. Check accessibility basics
        analysis['accessibility'] = await page.evaluate('''() => {
            const headingStructure = [];
            let lastLevel = 0;
            let skipCount = 0;

            document.querySelectorAll('h1, h2, h3, h4, h5, h6').forEach(h => {
                const level = parseInt(h.tagName.substring(1));
                if (level > lastLevel + 1) skipCount++;
                lastLevel = level;
                headingStructure.push(level);
            });

            return {
                hasLangAttribute: !!document.documentElement.lang,
                headingSkips: skipCount,
                altTextCoverage: Array.from(document.querySelectorAll('img')).filter(img => img.alt).length +
                                '/' + document.querySelectorAll('img').length,
                ariaLabels: document.querySelectorAll('[aria-label]').length
            };
        }''')

        # 9. Check mobile-friendliness
        await page.set_viewport_size({'width': 375, 'height': 667})  # iPhone size
        analysis['mobile'] = await page.evaluate('''() => {
            const body = document.body;
            const isScrollable = body.scrollWidth > window.innerWidth;
            const fontSize = window.getComputedStyle(body).fontSize;

            return {
                horizontalScroll: isScrollable,
                mobileFontSize: fontSize,
                viewportWidth: window.innerWidth
            };
        }''')

        await browser.close()

        return analysis

async def generate_critique(analysis):
    """Generate a critique based on the analysis"""

    critique = {
        'strengths': [],
        'weaknesses': [],
        'recommendations': []
    }

    # Typography analysis
    body_font_size = float(analysis['typography']['bodyFontSize'].replace('px', ''))
    if body_font_size >= 14:
        critique['strengths'].append(f"Good base font size ({analysis['typography']['bodyFontSize']})")
    else:
        critique['weaknesses'].append(f"Font size too small ({analysis['typography']['bodyFontSize']})")
        critique['recommendations'].append("Increase base font size to at least 16px for better readability")

    # Image analysis
    if analysis['images']['count'] > 0:
        imgs_without_alt = sum(1 for img in analysis['images']['sizes'] if not img['hasAlt'])
        if imgs_without_alt > 0:
            critique['weaknesses'].append(f"{imgs_without_alt} images lack alt text")
            critique['recommendations'].append("Add descriptive alt text to all images for accessibility")

    # Navigation
    if analysis['navigation']['hasTOC']:
        critique['strengths'].append("Has table of contents for navigation")
    else:
        critique['weaknesses'].append("No table of contents found")
        critique['recommendations'].append("Add a floating TOC for better navigation")

    # Interactivity
    if analysis['interactivity']['buttonCount'] == 0 and analysis['interactivity']['collapsibleCount'] == 0:
        critique['weaknesses'].append("No interactive elements found")
        critique['recommendations'].append("Add collapsible sections, tabs, or interactive charts")

    # Mobile
    if analysis['mobile']['horizontalScroll']:
        critique['weaknesses'].append("Content requires horizontal scrolling on mobile")
        critique['recommendations'].append("Ensure all content is responsive")

    # Accessibility
    if not analysis['accessibility']['hasLangAttribute']:
        critique['weaknesses'].append("Missing language attribute")
        critique['recommendations'].append("Add lang='en' to html element")

    if analysis['accessibility']['headingSkips'] > 0:
        critique['weaknesses'].append(f"Heading hierarchy has {analysis['accessibility']['headingSkips']} skip(s)")
        critique['recommendations'].append("Fix heading hierarchy for better screen reader navigation")

    return critique

async def main():
    """Main function to analyze HTML files"""

    html_files = [
        'output/manuscript/MANUSCRIPT_IMPROVED.html',
        'output/manuscript/COLLABORATOR_SUMMARY_IMPROVED.html',
        'output/manuscript/MANUSCRIPT_ULTRA_MODERN.html'
    ]

    results = {}

    for html_file in html_files:
        if Path(html_file).exists():
            print(f"\nAnalyzing {html_file}...")
            analysis = await analyze_html(html_file)
            critique = await generate_critique(analysis)

            results[html_file] = {
                'analysis': analysis,
                'critique': critique
            }

            # Print summary
            print(f"\n=== {html_file} ===")
            print(f"Strengths: {len(critique['strengths'])}")
            for s in critique['strengths']:
                print(f"  ✓ {s}")

            print(f"\nWeaknesses: {len(critique['weaknesses'])}")
            for w in critique['weaknesses']:
                print(f"  ✗ {w}")

            print(f"\nRecommendations: {len(critique['recommendations'])}")
            for r in critique['recommendations']:
                print(f"  → {r}")

    # Save full results
    with open('html_analysis_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print("\nFull analysis saved to html_analysis_results.json")

if __name__ == "__main__":
    asyncio.run(main())