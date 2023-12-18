import os

def getfiles(folder, extension, get_prefix=True):
    relative_paths = []
    prefixes = []

    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(extension):
                relative_path = os.path.relpath(os.path.join(root, file), folder)
                relative_paths.append(relative_path)
                if get_prefix:
                    # Extract everything until the last /
                    prefix = relative_path.replace("\\", "_").rstrip("_interact-ipro-ptm.pep.xml")
                    prefix = prefix.rstrip("_interact-prob") # handle cases where the file naming is different

                    # Strip the .zip extension if present
                    prefix, _ = os.path.splitext(prefix)

                    prefixes.append(prefix)

    return (relative_paths, prefixes) if get_prefix else relative_paths


