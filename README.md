# Puzzle with the Words "100 Days"

This program aims to use multiple images as a background while highlighting large images to form the word "100." 

**Note:** The program requires two sets of images:

1. **Background Image Set**: Contains multiple images to be used for generating the collage background.
2. **"100" Word Image Set**: Contains highlighted images to be used to form the word "100."

## Detailed Steps and Methods

### 1. Image Reading and Aspect Ratio Filtering

- **Getting Image Paths**: 
  - Use `os.listdir(image_dir)` to list all files in the specified directory and `os.path.join()` to create complete image paths.
  
- **Image Format Filtering**: 
  - Use list comprehension to filter supported image formats (e.g., PNG, JPG, GIF) to ensure only these file types are processed.

- **Aspect Ratio Calculation**: 
  - Open images using `Image.open(image_path)`, retrieve the width and height with `img.size`, and calculate the aspect ratio (width/height). 
  - Store each image's aspect ratio in the `aspect_ratios` list.

- **Filtering Conditions**: 
  - Filter images with `if aspect_ratio > 1.4:` to ensure they meet the desired criteria for the collage. 
  - Additionally, images with an aspect ratio less than 0.4 are removed by:
    - Using `shutil.move()` to move these images to a designated removal folder.

### 2. Generating the Background Image

- **Image Processing**: 
  - Use `img.resize((tile_size, tile_size))` to uniformly scale images for consistency.
  - Crop images with `img.crop()`:
    - For wider images, crop the center using `(left, 0, right, height)`.
    - For taller images, crop using `(0, top, width, bottom)`.

- **Creating the Background Canvas**: 
  - Use `Image.new('RGB', (tile_size * grid_size, tile_size * grid_size))` to create a new background image.

- **Collaging Images**: 
  - Use nested loops to iterate over the background canvas and apply `background.paste(img, (j * tile_size, i * tile_size))` to paste the images.

- **Opacity Adjustment**: 
  - Adjust the background image's opacity with `background.point(lambda p: p * 0.5)`, softening the background.

### 3. Adding Highlighted Text

- **Loading the Background Image**: 
  - Use `Image.open(background_path)` to load the generated background.

- **Processing the Text Image**: 
  - Read the images for the text and resize them with `img.resize((tile_size, tile_size))`.
  - Calculate the center position using `center_x` and `center_y` to ensure the text is centered.
  - Paste the text images onto the background with `background.paste(img, (paste_x, paste_y))`.

- **Saving and Displaying**: 
  - Save the final image using `background.save(output_image_path)` and display it with `background.show()`.

With these steps, the program can efficiently create collage images with highlighted text, suitable for memorial day composition and art creation.
