import dataclasses
import pytest

from clade_clustering.layout import Layout
from clade_clustering.paths import PROJECT_ROOT

# Only the properties which depends on the data version
VERSIONED_PROPERTIES = [
    "processed", "results", "manifests", 
    "sample_metadata_clean", "parsnp_vcf",
    "parsnp_matrices_dir", "mapping_panel_matrices_dir"
]


def test_v4_paths(tmp_path):
    layout = Layout("v4", root=tmp_path)

    assert layout.data == tmp_path / "data"

    assert layout.processed == tmp_path / "data" / "processed" / "v4"

    assert layout.results == tmp_path / "data" / "results" / "v4"

    assert layout.manifests == tmp_path / "manifests" / "v4"

    assert layout.sample_metadata_clean == tmp_path / "data" / "processed" / "v4" / "sample_metadata_clean.parquet"

    assert layout.parsnp_vcf == tmp_path / "data" / "results" / "v4" / "parsnp" / "parsnp.vcf"

    assert layout.parsnp_matrices_dir == tmp_path / "data" / "results" / "v4" / "matrices" / "parsnp"

    assert layout.mapping_panel_matrices_dir == tmp_path / "data" / "results" / "v4" / "matrices" / "mapping_panel"


def test_version_changes_versioned_paths(tmp_path):
    v4 = Layout("v4", root=tmp_path)
    v5 = Layout("v5", root=tmp_path)

    for name in VERSIONED_PROPERTIES:
        assert "v4" in getattr(v4, name).parts
        assert "v5" in getattr(v5, name).parts


def test_default_root_is_project_root():
    assert Layout("v4").root == PROJECT_ROOT


def test_layout_is_immutable(tmp_path):
    layout = Layout("v4", root=tmp_path)

    with pytest.raises(dataclasses.FrozenInstanceError):
        layout.version = "v5"

    # The frozen=True blocks attributions, before Python checks @property
    with pytest.raises(dataclasses.FrozenInstanceError):
        layout.processed = "x"


def test_layout_does_not_touch_disk(tmp_path):
    layout = Layout("v4", root=tmp_path)

    for name in ["data", *VERSIONED_PROPERTIES]:
        getattr(layout, name)

    assert list(tmp_path.iterdir()) == []