from pathlib import Path
import numpy as np


class GravityField:

    def __init__(self, source_file: Path) -> None:

        self.coeffs = np.loadtxt(source_file, skiprows=5).T

        return None


class General:

    def __init__(self, source_file: Path) -> None:

        header = False
        headers = {}
        count = 1
        current_header = {}
        current_data = []
        current_type = ""
        current = {}
        with source_file.open() as source:

            for line in source:
                if line.startswith("#"):
                    header = True
                    content = line.split(" ")[1:]
                    if content[0] == "name:":
                        current_header["id"] = content[1][:-1]
                    elif content[0] == "rows:":
                        current_header["rows"] = int(content[1])
                    elif content[0] == "columns:":
                        current_header["columns"] = int(content[1])
                    elif content[0] == "ndims:":
                        next(source)
                    elif content[0] == "type:":
                        current_header["type"] = content[1][:-1]
                else:
                    if header:
                        current_header["start"] = count + 1
                        header = False
                        current["header"] = current_header
                        current_header = {}
                    else:
                        if current["header"]["type"] in ["matrix", "vector", "scalar"]:
                            current_data.append(
                                np.array(line.split(" "), dtype=np.float64)
                            )

                count += 1

        for key, val in headers.items():
            print(key, val)

        return None
