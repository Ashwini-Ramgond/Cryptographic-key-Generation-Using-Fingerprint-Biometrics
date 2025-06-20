% --- Step 1: Fingerprint Minutiae Extraction ---
% Read Input Image
input_image = imread('101_1.tif');

% Convert to Grayscale (if not already)
if size(input_image, 3) == 3
    input_image = rgb2gray(input_image);
end

% Binarize the Image
binary_image = imbinarize(input_image);

% Thinning
thin_image = bwmorph(~binary_image, 'thin', Inf); % Use inverted binary for thinning
thin_image = ~thin_image; % Re-invert after thinning

% Initialize Variables for Minutiae Extraction
s = size(thin_image);
N = 3; % Window size
n = (N-1)/2;
temp = padarray(thin_image, [n, n], 0);
bifurcation = zeros(size(temp));
ridge = zeros(size(temp));

% Minutiae Extraction
for x = (n+1):(s(1)+n-1)
    for y = (n+1):(s(2)+n-1)
        mat = temp(x-n:x+n, y-n:y+n);
        if (mat(2, 2) == 0)
            ridge(x, y) = sum(sum(~mat)) == 2;
            bifurcation(x, y) = sum(sum(~mat)) == 4;
        end
    end
end

% Ridge Ending Detection
[ridge_x, ridge_y] = find(ridge == 1);

% Bifurcation Detection
[bifurcation_x, bifurcation_y] = find(bifurcation == 1);

% Combine Ridge Endings and Bifurcations
ridge_minutiae = [ridge_x, ridge_y, ones(size(ridge_x))]; % Ridge endings (type = 1)
bifurcation_minutiae = [bifurcation_x, bifurcation_y, zeros(size(bifurcation_x))]; % Bifurcations (type = 0)
all_minutiae = [ridge_minutiae; bifurcation_minutiae]; % Combine both

% Normalize Minutiae Coordinates
all_minutiae(:, 1) = all_minutiae(:, 1) / size(input_image, 1); % Normalize x
all_minutiae(:, 2) = all_minutiae(:, 2) / size(input_image, 2); % Normalize y

% Add Orientation as a Discriminative Feature
all_minutiae(:, 3) = all_minutiae(:, 3) + 0.5 * rand(size(all_minutiae(:, 3))); % Add small randomness

% --- Step 2: Bit String Generation ---
% Calculate Pairwise Distances Between Minutiae Points
distance_matrix = pdist2(all_minutiae(:, 1:3), all_minutiae(:, 1:3)); % Include type in distances

% Extract Unique Distances (Upper Triangle of Distance Matrix)
unique_distances = unique(triu(distance_matrix, 1));
unique_distances = unique_distances(unique_distances > 0.001); % Remove noise

% Sort Distances
sorted_distances = sort(unique_distances);

% Scale Distances to Increase Magnitude
scaling_factor = 1000; % Adjust this factor as needed
scaled_distances = round(sorted_distances * scaling_factor);

% Convert Distances to Binary and Then to Gray Code
combined_bit_string = '';
for i = 1:min(length(scaled_distances), 128) % Limit to first 128 distances
    % Convert Scaled Distance to Binary (10-bit representation)
    binary_distance = dec2bin(scaled_distances(i), 10);

    % Convert Binary to Gray Code
    gray_code = binary_distance;
    gray_code(1) = binary_distance(1);
    for j = 2:length(binary_distance)
        gray_code(j) = num2str(xor(str2double(binary_distance(j-1)), str2double(binary_distance(j))));
    end

    % Append to Combined Bit String
    combined_bit_string = [combined_bit_string, gray_code];
end

% Ensure the Bit String is of Fixed Length (1024 bits)
if length(combined_bit_string) > 1024
    combined_bit_string = combined_bit_string(1:1024); % Trim to 1024 bits
else
    combined_bit_string = pad(combined_bit_string, 1024, 'right', '0'); % Pad with zeros if shorter
end

% --- Step 3: Convert to Hexadecimal ---
% Convert the Bit String into Hexadecimal
hex_string = '';
for i = 1:4:length(combined_bit_string)
    % Take 4 bits at a time (since each hex digit corresponds to 4 bits)
    binary_group = combined_bit_string(i:min(i+3, length(combined_bit_string)));

    % Convert binary group to decimal and then to hexadecimal
    decimal_value = bin2dec(binary_group);
    hex_value = dec2hex(decimal_value);

    % Append the hex value to the final hex string
    hex_string = [hex_string, hex_value];
end

% Display the Generated Bit String
disp('Generated 1024-Bit String:');
disp(combined_bit_string);

% Display the Hexadecimal String
disp('Hexadecimal Representation of the Bit String:');
disp(hex_string);

% --- Visualize the Results ---
% Display the Images
figure;

% Display Input Image
subplot(1, 4, 1); imshow(input_image); title('Input Image');

% Display Binarized Image
subplot(1, 4, 2); imshow(binary_image); title('Binarized Image');

% Display Thinned Image
subplot(1, 4, 3); imshow(thin_image); title('Thinned Image');

% Create a white background
outImg = ones(size(thin_image, 1), size(thin_image, 2), 3) * 255; % White background (RGB)

% Highlight Ridge Endings (Red)
for i = 1:length(ridge_x)
    outImg(max(1, ridge_x(i)-2):min(size(outImg, 1), ridge_x(i)+2), ...
           max(1, ridge_y(i)-2):min(size(outImg, 2), ridge_y(i)+2), 1) = 255; % Red channel
    outImg(max(1, ridge_x(i)-2):min(size(outImg, 1), ridge_x(i)+2), ...
           max(1, ridge_y(i)-2):min(size(outImg, 2), ridge_y(i)+2), 2:3) = 0; % Remove green and blue
end

% Highlight Bifurcations (Blue)
for i = 1:length(bifurcation_x)
    outImg(max(1, bifurcation_x(i)-2):min(size(outImg, 1), bifurcation_x(i)+2), ...
           max(1, bifurcation_y(i)-2):min(size(outImg, 2), bifurcation_y(i)+2), 3) = 255; % Blue channel
    outImg(max(1, bifurcation_x(i)-2):min(size(outImg, 1), bifurcation_x(i)+2), ...
           max(1, bifurcation_y(i)-2):min(size(outImg, 2), bifurcation_y(i)+2), 1:2) = 0; % Remove red and green
end

% Display Final Image with Minutiae Highlighted
subplot(1, 4, 4); imshow(uint8(outImg)); title('Minutiae Highlighted');
% Extract Two Integers from the Bit String
integer1 = bin2dec(combined_bit_string(1:16)); % Use the first 16 bits
integer2 = bin2dec(combined_bit_string(17:32)); % Use the next 16 bits
% Ensure Both Integers are Prime
integer1 = max(2, integer1); % Ensure >= 2
integer2 = max(2, integer2); % Ensure >= 2
while ~isprime(integer1)
    integer1 = integer1 + 1;
end
while ~isprime(integer2)
    integer2 = integer2 + 1;
end


% Display the Two Integers
fprintf('Extracted Primes: p = %d, q = %d\n', integer1, integer2);

% Ensure the Primes are Set for RSA
p = integer1;  % Use the extracted prime integer1
q = integer2;  % Use the extracted prime integer2


disp('RSA algorithm');

n = p * q; % Calculate n
fprintf('\nn=%d', n);

phi = (p - 1) * (q - 1); % Calculate phi(n)
fprintf('\nphi(%d) is %d', n, phi);

val = 0;
cd = 0;

% Find e such that gcd(e, phi) = 1
while(cd ~= 1 || val == 0)
    n1 = randi(n, 1, 1);
    e = randi([2 n1], 1, 1);
    val = isprime(e);
    cd = gcd(e, phi);
end
fprintf('\ne = %d', e); % Displaying 'e'

val1 = 0;
d = 0;

% Find d such that d * e ≡ 1 (mod phi(n))
while(val1 ~= 1)
    d = d + 1;
    val1 = mod(d * e, phi);
end
fprintf('\nd=%d', d);
fprintf('\nPublic key is (%d,%d)', e, n);
fprintf('\nPrivate key is (%d,%d)', d, n);

% Input message and encrypt/decrypt
m = input('\nEnter the message: ', 's');
m1 = double(m); % ASCII equivalent of message
disp('ASCII equivalent of message ');
disp(m1);

over = length(m1);
o = 1;
c1 = zeros(1, over); % Initialize encryption vector
nm1 = char(zeros(1, over)); % Initialize decryption vector

% Encrypt message using RSA
while(o <= over)
    m = m1(o);
    diff = 0;
    if(m > n)
        diff = m - n + 1;
    end
    m = m - diff;

    c = 1;
    qm = dec2bin(e);
    len = length(qm);
    xz = 1;
    while(xz <= len)
        if(qm(xz) == '1')
            c = mod(mod((c^2), n) * m, n);
        elseif(qm(xz) == '0')
            c = (mod(c^2, n));
        end
        xz = xz + 1;
    end
    c1(o) = c;

    % Decrypt using d
    nm = 1;
    qm1 = dec2bin(d);
    len1 = length(qm1);
    xy = 1;
    while(xy <= len1)
        if(qm1(xy) == '1')
            nm = mod(mod((nm^2), n) * c, n);
        elseif(qm1(xy) == '0')
            nm = (mod(nm^2, n));
        end
        xy = xy + 1;
    end
    nm = nm + diff;
    nm1(o) = char(nm);
    o = o + 1;
end

% Display the encrypted message
o = 1;
fprintf('\nThe encrypted message is \n');
while(o <= over)
    fprintf('\t%d', c1(o));
    o = o + 1;
end

% Display the decrypted message in ASCII
o = 1;
fprintf('\nThe decrypted message in ASCII is \n');
while(o <= over)
    fprintf('\t%d', nm1(o));
    o = o + 1;
end

% Display the decrypted message
fprintf('\nThe decrypted message is: ');
disp(nm1);
fprintf('\n');
