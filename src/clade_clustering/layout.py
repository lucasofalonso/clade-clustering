from dataclasses import dataclass
from pathlib import Path

from clade_clustering.paths import PROJECT_ROOT


@dataclass(frozen=True)
class Layout:
    """
    Versioned directory layout and paths of the artifacts shared between pipeline stages.

    This class only resolves paths; it does not read, write or create anything on disk.

    Attributes:
        version: Dataset version processed by the pipeline (e.g. "v4").
        root: Project root directory.
    """
    version: str
    root: Path = PROJECT_ROOT

    # Main data directories
    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def processed(self) -> Path:
        return self.data / "processed" / self.version

    @property
    def results(self) -> Path:
        return self.data / "results" / self.version

    @property
    def manifests(self) -> Path:
        return self.root / "manifests" / self.version

    # Contracts between stages (producer -> consumer)
    @property
    def sample_metadata_clean(self) -> Path: # 01 -> 03
        return self.processed / "sample_metadata_clean.parquet"

    @property
    def parsnp_vcf(self) -> Path: # 02 -> 03
        return self.results / "parsnp" / "parsnp.vcf"

    @property
    def parsnp_matrices_dir(self) -> Path: # 03 -> 04
        return self.results / "matrices" / "parsnp" 

    @property
    def mapping_panel_matrices_dir(self) -> Path: # 03 -> 05
        return self.results / "matrices" / "mapping_panel"