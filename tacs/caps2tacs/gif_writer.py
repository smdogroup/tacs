"""
Written by Sean Engelstad, GT SMDO Lab, 2022-2023
"""
__all__ = ["GifWriter"]
import imageio, os


class GifWriter:
    """
    module to write gifs from a set of pngs
    """

    def __init__(self, duration: int = 20):
        # duration of each frame
        self._duration = duration

    def __call__(self, gif_filename: str, path: str):
        """
        call on current path to create gif of given filename
        """
        gif_filepath = os.path.join(path, gif_filename)
        with imageio.get_writer(
            gif_filepath, mode="I", duration=self._duration
        ) as writer:
            path_dir = os.listdir(path)
            path_dir = sorted(path_dir)
            for image_file in path_dir:
                print(image_file)
                if ".png" in image_file:
                    image = imageio.imread(image_file)
                    writer.append_data(image)


# example of how to use the GifWriter
# if __name__ == "__main__":
#    my_writer = GifWriter(frames_per_second=4)
#    my_writer("sizing1.gif", os.getcwd())
