from pathlib import Path
import sys
import tempfile


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import glossary_symbols as gs
import math_symbol_parser as msp
import paperlint_actions as actions
import tex_structure as ts


def test_preferred_terminology_pass_texts_do_not_change():
    pass_texts = [
        "The exterior region connects to the interior region.",
        "Exterior regions interact with damping plate hardware.",
    ]

    for text in pass_texts:
        updated, span = actions.preferred_terminology_line(text)
        assert updated == text
        assert span is None


def test_preferred_terminology_fail_texts_are_rewritten():
    fail_to_expected = {
        "The external region connects to the internal region.": "The exterior region connects to the interior region.",
        "These external regions interact with reaction plate hardware.": "These exterior regions interact with damping plate hardware.",
        "External region, internal region, reaction plate.": "Exterior region, interior region, damping plate.",
    }

    for original, expected in fail_to_expected.items():
        updated, span = actions.preferred_terminology_line(original)
        assert updated == expected
        assert span == (0, len(original))


def test_symbol_to_glossary_label_prime_and_star_cases():
    assert gs.symbol_to_glossary_label("p^{'}_F") == "sym-p-prime-F"
    assert gs.symbol_to_glossary_label(r"p_\prime_F") == "sym-p-prime-F"
    assert gs.symbol_to_glossary_label(r"a_{\mu^*\nu^*}") == "sym-a-mu-star-nu-star"
    assert gs.symbol_to_glossary_label(r"b_{\mu^*}") == "sym-b-mu-star"
    assert gs.symbol_to_glossary_label(r"\Gamma_{\mu^*}") == "sym-Gamma-mu-star"


def test_glossary_equivalence_collapses_text_and_braced_subscripts():
    assert gs._normalize_glossary_equivalence_symbol(r"F_\text{max}") == gs._normalize_glossary_equivalence_symbol(r"F_{max}")


def test_glossary_equivalence_prefers_styled_symbol():
    assert gs._prefer_glossary_symbol(r"F_{max}", r"F_\text{max}") == r"F_\text{max}"


def test_extract_math_symbols_preserves_braced_text_subscript():
    symbols = msp.extract_math_symbols(r"h_{1,\text{stiff},d}")
    assert r"h_{1,\text{stiff},d}" in symbols


def test_extract_math_symbols_preserves_braces_for_multi_char_subscript():
    symbols = msp.extract_math_symbols(r"a_{in}")
    assert r"a_{in}" in symbols


def test_missing_word_style_masks_gls_and_rewrites_plain_word():
    prose_style_for_word, styled_word_pattern = actions.collect_word_style_data(
        "Reference: \\textit{max} and max and max."
    )
    line = r"\gls{sym-max} and max should be styled."

    updated, span = actions.missing_word_style_line(
        line, prose_style_for_word, styled_word_pattern
    )

    assert r"\gls{sym-max}" in updated
    assert r"\textit{max}" in updated
    assert span == (0, len(line))


def test_missing_word_style_pass_text_already_styled():
    prose_style_for_word, styled_word_pattern = actions.collect_word_style_data(
        "Reference: \\textit{max} and max and max."
    )
    line = r"\gls{sym-max} and \textit{max} should stay styled."

    updated, span = actions.missing_word_style_line(
        line, prose_style_for_word, styled_word_pattern
    )

    assert updated == line
    assert span is None


def test_parse_glossary_entries_preserves_nested_gls_symbol_bodies():
    line1 = r"\glsxtrnewsymbol[description={},sort=sym-Lambda-gls-sym-k-index]{sym-vec-Lambda-gls-sym-k-index}{\ensuremath{\vec{\Lambda}_\gls{sym-k-index}}}"
    line2 = r"\glsxtrnewsymbol[description={}]{sym-d-gls-sym-mu-star-gls-sym-nu-star}{\ensuremath{d_{\gls{sym-mu}^*,\gls{sym-nu}^*}}}"

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "symbols.tex"
        path.write_text(line1 + "\n" + line2 + "\n")
        entries = msp.parse_glossary_entries(path)

    symbols = {symbol for _, symbol, _ in entries}
    assert r"\vec{\Lambda}_\gls{sym-k-index}" in symbols
    assert r"d_{\gls{sym-mu}^*,\gls{sym-nu}^*}" in symbols


def test_collect_tex_files_finds_nested_includes_with_root_fallback():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        main = root / "main.tex"
        section = root / "sections" / "sec1.tex"
        grand = root / "appendix-note.tex"

        section.parent.mkdir(parents=True)
        main.write_text("\\input{sections/sec1}\n")
        section.write_text("\\input{appendix-note}\n")
        grand.write_text("Nested include target.\n")

        files = ts.collect_tex_files(str(main))

    file_names = {Path(item).name for item in files}
    assert "main.tex" in file_names
    assert "sec1.tex" in file_names
    assert "appendix-note.tex" in file_names
