Cryptographic-key-Generation-Using-Fingerprint-Biometrics

This project focuses on generating a secure cryptographic key using biometric fingerprint data. Implemented using MATLAB, this method ensures enhanced security by integrating unique biometric features for RSA encryption.
Project Overview
Traditional cryptographic systems rely on password-based or randomly generated keys. However, these can be stolen, guessed, or lost. This project proposes a secure alternative: generating encryption keys using biometric fingerprints, which are unique to individuals.

Objectives

Extract unique features (minutiae) from fingerprint images

Convert extracted data into a binary sequence

Generate a 128-bit cryptographic key

Use this key for RSA encryption and decryption

Simulate the entire process using MATLAB

Tools & Technologies

MATLAB (Fingerprint image processing, key generation)

MATLAB Image Processing Toolbox

RSA Algorithm (for encryption/decryption)

FPGA (Spartan-6) (for hardware implementation — in progress)

Methodology

Fingerprint Acquisition: Load a fingerprint image (BMP format)

Preprocessing: Binarization, noise removal, thinning

Minutiae Extraction: Identify ridge endings and bifurcations

Key Generation: Convert minutiae to binary and derive a 128-bit key

Encryption/Decryption: Use RSA algorithm with generated key

(Planned): Convert MATLAB code to Verilog HDL for FPGA implementation
 

User input: fingerprint image file

Output: binary key, encrypted and decrypted messages

Results

Successfully generated 1024-bit biometric keys

RSA encryption and decryption verified in both MATLAB and FPGA simulation

Keys dynamically generated from fingerprints—no static storage required

Verified simulation outputs using Xilinx ISE tool


Benefits

Eliminates need for remembering passwords

Provides individual-specific encryption

Resistant to key theft or duplication

Future Scope

Real-time fingerprint capture using sensors

Implementation on FPGA using Verilog HDL

Integration with secure digital systems (e.g., voting, banking)

Extendable to other biometric traits (iris, face)


IEEE papers on Biometric Key Generation

MATLAB documentation on Image Processing

RSA Algorithm Theory and Applications

