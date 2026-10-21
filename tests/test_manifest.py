from pathlib import Path

from arrayevidence.manifest import read_assays


def test_manifest_sections_and_coordinates(tmp_path: Path):
    manifest = tmp_path / "synthetic_manifest.csv"

    manifest.write_text(
        "[Heading]\n"
        "Descriptor,Value\n"
        "[Assay]\n"
        "IlmnID,Name,GenomeBuild,Chr,MapInfo,SNP,SourceSeq\n"
        "probe_A,rs100,38,1,1000,[A/G],AAA[A/G]CCC\n"
        "probe_B,rs101,38,1,1000,[C/T],AAA[C/T]CCC\n"
        "probe_C,unmapped,38,0,0,[A/C],AAA[A/C]CCC\n"
        "[Controls]\n"
        "Control,Value\n"
        "control_A,123\n",
        encoding="utf-8",
    )

    assays = list(read_assays(manifest))

    assert len(assays) == 3
    assert assays[0].coordinate == ("1", 1000)
    assert assays[1].coordinate == ("1", 1000)
    assert assays[2].coordinate is None
    assert assays[0].name == "rs100"
