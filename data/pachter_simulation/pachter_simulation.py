#!/usr/bin/env python

# Made by Paul Kiessling pakiessling@ukaachen.de


import os
import argparse
import tempfile
import shutil
from pypdl import Downloader

LINKS = {
    "https://zenodo.org/records/13944111/files/concordex_sim.zip":"fb7c79fd9cec2c79e3b74fb50be40ff4"
}


def download_links(links, temp_dir):
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:47.0) Gecko/20100101 Firefox/47.0"
    }
    dl = Downloader(headers=headers)
    for link, checksum in links.items():
        print(f"Downloading {link}")
        file = dl.start(
            url=link,
            file_path=temp_dir,
            segments=10,
            display=True,
            multithread=True,
            block=True,
            retries=3,
        )
        if not file.validate_hash(checksum, "md5"):
            raise ValueError(f"File {file} is corrupted")




def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="Convert Visium HD data to Spacehack format."
    )

    # Add arguments for output folder
    parser.add_argument(
        "-o", "--out_dir", help="Output directory to write files to.", required=True
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Download and process
    with tempfile.TemporaryDirectory() as temp_dir:  #
        download_links(LINKS, temp_dir)
        for file in os.listdir(temp_dir):
            if file.endswith(".tar.gz") or file.endswith(".zip"):
                shutil.unpack_archive(os.path.join(temp_dir, file), args.out_dir)
           


if __name__ == "__main__":
    main()
