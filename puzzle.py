import os
from PIL import Image
import matplotlib.pyplot as plt
import shutil
import random


# 绘制长宽比的条形统计图
image_dir = "/home/mahc/data/picture"
image_paths = [os.path.join(image_dir, img) for img in os.listdir(image_dir) if img.endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif'))]

# 计算每张图片的长宽比（宽/高）
aspect_ratios = []
for image_path in image_paths:
    try:
        with Image.open(image_path) as img:
            width, height = img.size
            aspect_ratios.append(width / height)
    except Exception as e:
        print(f"Error processing image {image_path}: {e}")

plt.figure(figsize=(10, 6))
counts, bins, patches = plt.hist(aspect_ratios, bins=30, color='skyblue', edgecolor='black')

for count, bin, patch in zip(counts, bins, patches):
    height = patch.get_height()
    if height > 0: 
        plt.text(patch.get_x() + patch.get_width() / 2, height, int(height), ha='center', va='bottom')

plt.title('Aspect Ratio Distribution of Images')
plt.xlabel('Aspect Ratio (Width / Height)')
plt.ylabel('Number of Images')
plt.grid(True)
plt.show()



#筛选照片，自己调阈值
image_dir = "/home/mahc/data/picture"
output_dir = "/home/mahc/data/filtered_images"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 筛选出长宽比大于1.4的图片
filtered_images = []
for image_path in image_paths:
    try:
        with Image.open(image_path) as img:
            width, height = img.size
            aspect_ratio = width / height
            if aspect_ratio > 1.4:
                filtered_images.append(image_path)
    except Exception as e:
        print(f"Error processing image {image_path}: {e}")

# 输出
for i, img_path in enumerate(filtered_images[:6]):
    try:
        with Image.open(img_path) as img:
            img.save(os.path.join(output_dir, f"filtered_image_{i+1}.jpg"))  
            plt.figure()
            plt.imshow(img)
            plt.title(f"Image path: {img_path}, Aspect ratio: {img.size[0] / img.size[1]:.2f}")
            plt.axis('off') 
            plt.show()
            
            print(f"Image path: {img_path}, Aspect ratio: {img.size[0] / img.size[1]}")
    except Exception as e:
        print(f"Error processing image {img_path}: {e}")



if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 筛选出长宽比小于0.4的图片
filtered_images = []
for image_path in image_paths:
    try:
        with Image.open(image_path) as img:
            width, height = img.size
            aspect_ratio = width / height
            if aspect_ratio < 0.4:
                filtered_images.append(image_path)
    except Exception as e:
        print(f"Error processing image {image_path}: {e}")

for i, img_path in enumerate(filtered_images[:6]):
    try:
        with Image.open(img_path) as img:
            img.save(os.path.join(output_dir, f"filtered_image_{i+1}.jpg"))  
            plt.figure()
            plt.imshow(img)
            plt.title(f"Image path: {img_path}, Aspect ratio: {img.size[0] / img.size[1]:.2f}")
            plt.axis('off')
            plt.show()
            
            print(f"Image path: {img_path}, Aspect ratio: {img.size[0] / img.size[1]}")
    except Exception as e:
        print(f"Error processing image {img_path}: {e}")

#移除小于0.4的
image_dir = "/home/mahc/data/picture"
remove_dir = "/home/mahc/data/removed_images"

if not os.path.exists(remove_dir):
    os.makedirs(remove_dir)

image_paths = [os.path.join(image_dir, img) for img in os.listdir(image_dir) if img.endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif'))]

for image_path in image_paths:
    try:
        with Image.open(image_path) as img:
            width, height = img.size
            aspect_ratio = width / height
            if aspect_ratio < 0.4:
                shutil.move(image_path, os.path.join(remove_dir, os.path.basename(image_path)))
                print(f"Moved image: {image_path} to {remove_dir}")
    except Exception as e:
        print(f"Error processing image {image_path}: {e}")


# 打印统计结果

image_dir = "/home/mahc/data/picture"
image_paths = [os.path.join(image_dir, img) for img in os.listdir(image_dir) if img.endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif'))]
total_images = len(image_paths)
portrait_images = 0  
landscape_images = 0  

for image_path in image_paths:
    try:
        with Image.open(image_path) as img:
            width, height = img.size
            if width > height:
                landscape_images += 1  
            else:
                portrait_images += 1 
    except Exception as e:
        print(f"Error processing image {image_path}: {e}")

print(f"Total images: {total_images}")
print(f"Portrait images (宽比长大): {portrait_images}")
print(f"Landscape images (长比宽大): {landscape_images}")



#背景图片
image_dir = "/home/mahc/data/picture"
output_image_path = "/home/mahc/data/background_200.jpg"

image_paths = [os.path.join(image_dir, img) for img in os.listdir(image_dir) if img.endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif'))]


tile_size = 1000   #像素
grid_size = 20     #背景长宽


random.shuffle(image_paths)


processed_images = []
for image_path in image_paths:
    try:
        with Image.open(image_path) as img:
            width, height = img.size
            if width > height:
                left = (width - height) // 2
                right = left + height
                img = img.crop((left, 0, right, height))
            elif height > width:
                top = (height - width) // 2
                bottom = top + width
                img = img.crop((0, top, width, bottom))

         
            img = img.resize((tile_size, tile_size))
            processed_images.append(img)
    except Exception as e:
        print(f"Error processing image {image_path}: {e}")


background = Image.new('RGB', (tile_size * grid_size, tile_size * grid_size))

for i in range(grid_size):
    for j in range(grid_size):
        if i * grid_size + j < len(processed_images):
            img = processed_images[i * grid_size + j]
            background.paste(img, (j * tile_size, i * tile_size))


background = background.point(lambda p: p * 0.5)

background.save(output_image_path)
background.show()



#加上100字样

Image.MAX_IMAGE_PIXELS = None  # 取消限制

image_dir = "/home/mahc/data/100/100"
background_path = "/home/mahc/data/background_200.jpg"
output_image_path = "/home/mahc/data/output_100.jpg"

image_paths = [os.path.join(image_dir, img) for img in os.listdir(image_dir) if img.endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif'))]

random.shuffle(image_paths)
image_paths = image_paths[:53]

background = Image.open(background_path)

tile_size = 1700  # 自己调
bg_width, bg_height = background.size

total_width = 9 * tile_size 
total_height = 9 * tile_size  

center_x = (bg_width - total_width) // 2
center_y = (bg_height - total_height) // 2

layout_1 = [(0, i) for i in range(9)]  


layout_0_1 = [(2, j) for j in range(9)] + [(3, 0), (3, 1), (3, 7), (3, 8)] + [(4, j) for j in range(9)]


layout_0_2 = [(6, j) for j in range(9)] + [(7, 0), (7, 1), (7, 7), (7, 8)] + [(8, j) for j in range(9)]


layout = layout_1 + layout_0_1 + layout_0_2


for idx, (x, y) in enumerate(layout):
    img = Image.open(image_paths[idx])
    
    width, height = img.size
    if width > height:
        img = img.crop(((width - height) // 2, 0, (width + height) // 2, height))
    else:
        img = img.crop((0, (height - width) // 2, width, (height + width) // 2))
    
    img = img.resize((tile_size, tile_size))

    paste_x = center_x + x * tile_size
    paste_y = center_y + y * tile_size
    
    background.paste(img, (paste_x, paste_y))

background.save(output_image_path)
background.show()
