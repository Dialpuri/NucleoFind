import xml.etree.ElementTree as ET
from pathlib import Path
import nucleofind as nf


def test_doc_version():
    current_file = Path(__file__)
    base_dir = current_file.parents[2]
    doc_path = base_dir / "docs" / "Writerside" / "writerside.cfg"

    assert doc_path.exists()
    xml = ET.parse(str(doc_path))
    root = xml.getroot()
    x = root.find('instance')
    doc_version = x.attrib['version']
    assert nf.__version__ == doc_version
