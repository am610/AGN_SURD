"""Retrieve the exact public inputs used by the adopted analysis."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.error
import urllib.request

HERE = Path(__file__).resolve().parent


def retrieve(url):
    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            return response.read()
    except (urllib.error.URLError, TimeoutError):
        if not shutil.which('curl'):
            raise
        return subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error',
                               '--max-time', '30', url], check=True, capture_output=True).stdout


def install(path, data, digest):
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f'Source hash mismatch: {path.name}')
    if path.exists():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'Existing file differs; preserved without overwriting: {path}')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=HERE.parent)
    args = parser.parse_args()
    config = json.loads((HERE / 'input_sources.json').read_text())
    for source in config['files']:
        path = args.root / source['path']
        if path.exists():
            data = path.read_bytes()
        else:
            data = retrieve(source['url'])
            if source['encoding'] == 'gzip':
                data = gzip.decompress(data)
        install(path, data, source['sha256'])
        print(f'Verified {source["path"]}', flush=True)
    folder = args.root / 'agn_surd_project/agn_data/ngc5548_agnwatch'
    with tarfile.open(fileobj=io.BytesIO((folder / 'hb_profiles.tar.gz').read_bytes()), mode='r:gz') as archive:
        members = {m.name: m for m in archive.getmembers() if m.isfile()}
        if set(members) != set(config['spectra_sha256']):
            raise ValueError('Unexpected spectrum archive contents')
        for name, digest in config['spectra_sha256'].items():
            if Path(name).name != name:
                raise ValueError('Invalid spectrum name')
            with archive.extractfile(members[name]) as stream:
                install(folder / 'hb_profiles_extracted' / name, stream.read(), digest)
    print(f'Verified {len(members)} spectra and the pinned SURD implementation')


if __name__ == '__main__':
    main()
