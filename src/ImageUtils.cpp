//
// Created by Sofia Iannicelli on 2/18/25.
//

#include "ImageUtils.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../include/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../include/stb_image_write.h"

unsigned char * ImageUtils::load_img(const char *filename, int *width, int *height,
                                     int *channels_in_file, int desired_channels) {
    return stbi_load(filename, width, height, channels_in_file, desired_channels);
}

int ImageUtils::write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes) {
    return stbi_write_png(filename, w, h, comp, data, stride_in_bytes);
}
