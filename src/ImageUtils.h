//
// Created by Sofia Iannicelli on 2/18/25.
//

#ifndef RAYTRACER_2_IMAGEUTILS_H
#define RAYTRACER_2_IMAGEUTILS_H

/**
 * Wrapper for stbi_image loading and writing
 */
class ImageUtils {
public:
    /**
     * Load image intensity values into array
     * @param filename file where image is stored
     * @param[out] width returns pixel width
     * @param[out] height returns pixel height
     * @param[out] channels_in_file returns total channels
     * @param desired_channels how many channels to extract
     * @return array of image values
     */
    static unsigned char *load_img(const char* filename, int *width, int *height,
                                   int *channels_in_file, int desired_channels);

    /**
     *
     * @param filename
     * @param w width in number of pixels
     * @param h height in number of pixels
     * @param comp # of components per pixel
     * @param data pointer to first byte of top-left-most pixel
     * @param stride_in_bytes distance in bytes from the first byte of
     * a row of pixels to the first byte of the next row of pixels
     * @return 0 on failure and non-0 on success
     */
    static int write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes);

};


#endif //RAYTRACER_2_IMAGEUTILS_H
