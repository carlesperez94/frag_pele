import os


def create_symlinks(src, dst):
    if not os.path.islink(dst):
        os.symlink(src, dst)
