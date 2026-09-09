"""Tests for the prompt registry: strict binding, gating, and provenance."""

import pytest

from panspatial.orchestrator.registry import PromptError, PromptRegistry


@pytest.fixture(scope="module")
def registry():
    return PromptRegistry("prompts")


def test_every_prompt_loads_and_the_system_anchor_exists(registry):
    assert "system.meta" in registry.ids()
    assert len(registry) >= 7
    assert "unit of replication is the PATIENT" in registry.system_prompt


def test_task_prompts_all_declare_their_placeholders(registry):
    """A placeholder not declared in frontmatter can never be bound, so it would ship raw."""
    for tid in registry.ids():
        template = registry.get(tid)
        undeclared = template.declared_placeholders() - set(template.variables)
        assert not undeclared, f"{tid} uses undeclared placeholders: {sorted(undeclared)}"


def test_missing_required_variable_is_refused(registry):
    with pytest.raises(PromptError, match="missing required variables"):
        registry.render("phase3.nk_spatial", {"adata_path": "x.h5ad"})


def test_undeclared_variable_is_refused(registry):
    with pytest.raises(PromptError, match="undeclared variables"):
        registry.render(
            "phase3.nk_spatial",
            {"adata_path": "x.h5ad", "platform": "Xenium", "typo_var": "1"},
        )


def test_defaults_fill_optional_variables(registry):
    rendered = registry.render("phase3.nk_spatial", {"adata_path": "x.h5ad", "platform": "Xenium"})
    assert "malignant, mCAF, iCAF, TLS, vasculature" in rendered.user
    assert "{{" not in rendered.user


def test_manuscript_prompts_are_gated_on_verified_results(registry):
    values = {
        "core_finding": "CAF-NK axis",
        "results_table": "see table",
        "n_tumor_types": "12",
        "limitations": "observational",
    }
    with pytest.raises(PromptError, match="gated"):
        registry.render("phase5.title_abstract", values)
    rendered = registry.render("phase5.title_abstract", values, results_verified=True)
    assert "CAF-NK axis" in rendered.user


def test_system_anchor_cannot_be_run_as_a_task(registry):
    with pytest.raises(PromptError, match="system anchor"):
        registry.render("system.meta", {})


def test_unknown_prompt_id_lists_the_available_ones(registry):
    with pytest.raises(PromptError, match="phase3.nk_spatial"):
        registry.get("phase3.nk_spatiaI")


def test_fingerprint_is_stable_and_input_sensitive(registry):
    a = registry.render("phase3.nk_spatial", {"adata_path": "x.h5ad", "platform": "Xenium"})
    b = registry.render("phase3.nk_spatial", {"adata_path": "x.h5ad", "platform": "Xenium"})
    c = registry.render("phase3.nk_spatial", {"adata_path": "y.h5ad", "platform": "Xenium"})
    d = registry.render("phase3.nk_spatial", {"adata_path": "x.h5ad", "platform": "Xenium"}, effort="low")
    assert a.fingerprint == b.fingerprint
    assert a.fingerprint != c.fingerprint
    assert a.fingerprint != d.fingerprint


def test_system_anchor_is_attached_to_every_render(registry):
    rendered = registry.render("phase6.audit", {"target_paths": "a.py", "source_bundle": "pass"})
    assert rendered.system == registry.system_prompt
    assert rendered.model == "claude-opus-5"
