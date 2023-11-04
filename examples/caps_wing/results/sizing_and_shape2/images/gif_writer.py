import imageio, os


class GifWriter:
    """
    module to write gifs from a set of pngs
    """

    def __init__(self, frames_per_second: int = 4):
        self._fps = frames_per_second

    def __call__(self, gif_filename: str, path: str):
        """
        call on current path to create gif of given filename
        """
        gif_filepath = os.path.join(path, gif_filename)
        with imageio.get_writer(gif_filepath, mode="I", fps=self._fps) as writer:
            path_dir = os.listdir(path)
            path_dir = sorted(path_dir)
            for image_file in path_dir:
                print(image_file)
                if ".png" in image_file:
                    image = imageio.imread(image_file)
                    writer.append_data(image)


my_writer = GifWriter(frames_per_second=20)
my_writer("size_and_shape2.gif", os.getcwd())
