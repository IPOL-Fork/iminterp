#! /bin/bash
# TV denoising example BASH shell script

# Echo shell commands
set -v

# Generate Gaussian noise with standard deviation 15 on "einstein.bmp"
# and save the result to "noisy.bmp".
./imnoise gaussian:15 einstein.bmp noisy.bmp

# Perform TV regularized denoising with the split Bregman algorithm on 
# "noisy.bmp" and save the result to "denoised.bmp".
./tvdenoise -n gaussian:15 noisy.bmp denoised.bmp

# Compare the original to the denoised image
./imdiff einstein.bmp denoised.bmp

