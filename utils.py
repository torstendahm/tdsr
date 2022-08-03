from functools import partial
import http.server
import shutil
import webbrowser
from pathlib import Path, PurePath
import socketserver
from tdsm.types import PathLike
import pprint
from typing import List


def delete_dir(f: Path):
    shutil.rmtree(f, ignore_errors=True)


def delete_file(f: Path):
    try:
        f.unlink(missing_ok=True)
    except TypeError:
        # missing_ok argument added in 3.8
        try:
            f.unlink()
        except FileNotFoundError:
            pass


def is_relative_to(p: PathLike, *other: List[PathLike]):
    try:
        p.relative_to(*other)
        return True
    except ValueError:
        return False


def serve_dir(directory: PathLike, port=8080):
    handler = partial(StaticFileHandler, base_dir=directory)
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.TCPServer(("", port), handler) as httpd:
        webbrowser.open(f"http://localhost:{port}")
        httpd.serve_forever()


class StaticFileHandler(http.server.BaseHTTPRequestHandler):
    def __init__(
        self, request: str, client_address: str, server: str, base_dir: PathLike = None
    ):
        self.base_dir = base_dir
        super().__init__(request, client_address, server)

    def send_error(self, error: Exception):
        if isinstance(error, FileNotFoundError):
            self.send_message(404, "not found")
        elif isinstance(error, PermissionError):
            self.send_message(403, "permission denied")
        else:
            self.send_message(500, str(error))

    def send_message(self, status: int, message: str):
        self.send_response(status)
        self.end_headers()
        self.wfile.write(message.encode())

    def send_dir(self, path: PathLike):
        contents = pprint.pformat([str(p) for p in path.iterdir()], indent=4)
        self.send_message(200, contents)

    def send_file(self, path: PathLike):
        with open(path, "rb") as f:
            data = f.read()
            self.send_response(200)
            self.end_headers()
            self.wfile.write(data)

    def do_GET(self):
        try:
            path = Path(self.path)
            path = self.base_dir / path.relative_to(path.anchor)
            if not is_relative_to(path, self.base_dir):
                self.send_message(403, "forbidden")
            elif (path / "index.html").is_file():
                self.send_file(path / "index.html")
            elif path.is_dir():
                self.send_dir(path)
            else:
                self.send_file(path)
        except Exception as e:
            self.send_error(e)
