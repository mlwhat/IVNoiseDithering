#define _CRT_SECURE_NO_WARNINGS

#include <windows.h>  // for bitmap headers.  Sorry non windows people!
#include <stdint.h>
#include <vector>
#include <random>
#include <atomic>
#include <thread>
#include <complex>
#include <array>

typedef uint8_t uint8;
typedef int64_t int64;

const float c_pi = 3.14159265359f;

// settings
const bool c_doDFT = true;

// globals 
FILE* g_logFile = nullptr;

//======================================================================================
inline float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
}

//======================================================================================
struct SColor
{
	SColor(uint8 _R = 0, uint8 _G = 0, uint8 _B = 0)
		: R(_R), G(_G), B(_B)
	{

	}

	inline void Set(uint8 _R, uint8 _G, uint8 _B)
	{
		R = _R;
		G = _G;
		B = _B;
	}

	uint8 B, G, R;
};

//======================================================================================
struct SImageData
{
	SImageData()
		: m_width(0)
		, m_height(0)
	{ }

	SColor& GetPixel(size_t x, size_t y)
	{
		return *(SColor*)&m_pixels[y * m_pitch + x * 3];
	}

	const SColor& GetPixel(size_t x, size_t y) const
	{
		return *(SColor*)&m_pixels[y * m_pitch + x * 3];
	}

	SColor& GetPixel(size_t index)
	{
		size_t y = index / m_width;
		size_t x = index % m_width;
		return GetPixel(x, y);
	}

	const SColor& GetPixel(size_t index) const
	{
		size_t y = index / m_width;
		size_t x = index % m_width;
		return GetPixel(x, y);
	}

	size_t m_width;
	size_t m_height;
	size_t m_pitch;
	std::vector<uint8> m_pixels;
};



//======================================================================================
struct SImageDataComplex
{
	SImageDataComplex()
		: m_width(0)
		, m_height(0)
	{ }

	size_t m_width;
	size_t m_height;
	std::vector<std::complex<float>> m_pixels;
};

//======================================================================================
std::complex<float> DFTPixel(const SImageData &srcImage, size_t K, size_t L)
{
	std::complex<float> ret(0.0f, 0.0f);

	for (size_t x = 0; x < srcImage.m_width; ++x)
	{
		for (size_t y = 0; y < srcImage.m_height; ++y)
		{
			// Get the pixel value (assuming greyscale) and convert it to [0,1] space
			const uint8 *src = &srcImage.m_pixels[(y * srcImage.m_pitch) + x * 3];
			float grey = float(src[0]) / 255.0f;

			// Add to the sum of the return value
			float v = float(K * x) / float(srcImage.m_width);
			v += float(L * y) / float(srcImage.m_height);
			ret += std::complex<float>(grey, 0.0f) * std::polar<float>(1.0f, -2.0f * c_pi * v);
		}
	}

	return ret;
}

//======================================================================================
void ImageDFT(const SImageData &srcImage, SImageDataComplex &destImage)
{
	// NOTE: this function assumes srcImage is greyscale, so works on only the red component of srcImage.
	// ImageToGrey() will convert an image to greyscale.

	// size the output dft data
	destImage.m_width = srcImage.m_width;
	destImage.m_height = srcImage.m_height;
	destImage.m_pixels.resize(destImage.m_width*destImage.m_height);

	size_t numThreads = std::thread::hardware_concurrency();
	//if (numThreads > 0)
		//numThreads = numThreads - 1;

	std::vector<std::thread> threads;
	threads.resize(numThreads);

	printf("Doing DFT with %zu threads...\n", numThreads);

	// calculate 2d dft (brute force, not using fast fourier transform) multithreadedly
	std::atomic<size_t> nextRow(0);
	for (std::thread& t : threads)
	{
		t = std::thread(
			[&]()
		{
			size_t row = nextRow.fetch_add(1);
			bool reportProgress = (row == 0);
			int lastPercent = -1;

			while (row < srcImage.m_height)
			{
				// calculate the DFT for every pixel / frequency in this row
				for (size_t x = 0; x < srcImage.m_width; ++x)
				{
					destImage.m_pixels[row * destImage.m_width + x] = DFTPixel(srcImage, x, row);
				}

				// report progress if we should
				if (reportProgress)
				{
					int percent = int(100.0f * float(row) / float(srcImage.m_height));
					if (lastPercent != percent)
					{
						lastPercent = percent;
						printf("            \rDFT: %i%%", lastPercent);
					}
				}

				// go to the next row
				row = nextRow.fetch_add(1);
			}
		}
		);
	}

	for (std::thread& t : threads)
		t.join();

	printf("\n");
}

//======================================================================================
void GetMagnitudeData(const SImageDataComplex& srcImage, SImageData& destImage)
{
	// size the output image
	destImage.m_width = srcImage.m_width;
	destImage.m_height = srcImage.m_height;
	destImage.m_pitch = 4 * ((srcImage.m_width * 24 + 31) / 32);
	destImage.m_pixels.resize(destImage.m_pitch*destImage.m_height);

	// get floating point magnitude data
	std::vector<float> magArray;
	magArray.resize(srcImage.m_width*srcImage.m_height);
	float maxmag = 0.0f;
	for (size_t x = 0; x < srcImage.m_width; ++x)
	{
		for (size_t y = 0; y < srcImage.m_height; ++y)
		{
			// Offset the information by half width & height in the positive direction.
			// This makes frequency 0 (DC) be at the image origin, like most diagrams show it.
			int k = (x + (int)srcImage.m_width / 2) % (int)srcImage.m_width;
			int l = (y + (int)srcImage.m_height / 2) % (int)srcImage.m_height;
			const std::complex<float> &src = srcImage.m_pixels[l*srcImage.m_width + k];

			float mag = std::abs(src);
			if (mag > maxmag)
				maxmag = mag;

			magArray[y*srcImage.m_width + x] = mag;
		}
	}
	if (maxmag == 0.0f)
		maxmag = 1.0f;

	const float c = 255.0f / log(1.0f + maxmag);

	// normalize the magnitude data and send it back in [0, 255]
	for (size_t x = 0; x < srcImage.m_width; ++x)
	{
		for (size_t y = 0; y < srcImage.m_height; ++y)
		{
			float src = c * log(1.0f + magArray[y*srcImage.m_width + x]);

			uint8 magu8 = uint8(src);

			uint8* dest = &destImage.m_pixels[y*destImage.m_pitch + x * 3];
			dest[0] = magu8;
			dest[1] = magu8;
			dest[2] = magu8;
		}
	}
}

//======================================================================================
bool ImageSave(const SImageData &image, const char *fileName)
{
	// open the file if we can
	FILE *file;
	file = fopen(fileName, "wb");
	if (!file) {
		printf("Could not save %s\n", fileName);
		return false;
	}

	// make the header info
	BITMAPFILEHEADER header;
	BITMAPINFOHEADER infoHeader;

	header.bfType = 0x4D42;
	header.bfReserved1 = 0;
	header.bfReserved2 = 0;
	header.bfOffBits = 54;

	infoHeader.biSize = 40;
	infoHeader.biWidth = (LONG)image.m_width;
	infoHeader.biHeight = (LONG)image.m_height;
	infoHeader.biPlanes = 1;
	infoHeader.biBitCount = 24;
	infoHeader.biCompression = 0;
	infoHeader.biSizeImage = (DWORD)image.m_pixels.size();
	infoHeader.biXPelsPerMeter = 0;
	infoHeader.biYPelsPerMeter = 0;
	infoHeader.biClrUsed = 0;
	infoHeader.biClrImportant = 0;

	header.bfSize = infoHeader.biSizeImage + header.bfOffBits;

	// write the data and close the file
	fwrite(&header, sizeof(header), 1, file);
	fwrite(&infoHeader, sizeof(infoHeader), 1, file);
	fwrite(&image.m_pixels[0], infoHeader.biSizeImage, 1, file);
	fclose(file);

	return true;
}

//======================================================================================
bool ImageLoad(const char *fileName, SImageData& imageData)
{
	// open the file if we can
	FILE *file;
	file = fopen(fileName, "rb");
	if (!file)
		return false;

	// read the headers if we can
	BITMAPFILEHEADER header;
	BITMAPINFOHEADER infoHeader;
	if (fread(&header, sizeof(header), 1, file) != 1 ||
		fread(&infoHeader, sizeof(infoHeader), 1, file) != 1 ||
		header.bfType != 0x4D42 || infoHeader.biBitCount != 24)
	{
		fclose(file);
		return false;
	}

	// read in our pixel data if we can. Note that it's in BGR order, and width is padded to the next power of 4
	imageData.m_pixels.resize(infoHeader.biSizeImage);
	fseek(file, header.bfOffBits, SEEK_SET);
	if (fread(&imageData.m_pixels[0], imageData.m_pixels.size(), 1, file) != 1)
	{
		fclose(file);
		return false;
	}

	imageData.m_width = infoHeader.biWidth;
	imageData.m_height = infoHeader.biHeight;
	imageData.m_pitch = 4 * ((imageData.m_width * 24 + 31) / 32);

	fclose(file);
	return true;
}

//======================================================================================
void ImageInit(SImageData& image, size_t width, size_t height)
{
	image.m_width = width;
	image.m_height = height;
	image.m_pitch = 4 * ((width * 24 + 31) / 32);
	image.m_pixels.resize(image.m_pitch * image.m_height);
	std::fill(image.m_pixels.begin(), image.m_pixels.end(), 255);
}

//======================================================================================
template <typename LAMBDA>
void ImageForEachPixel(SImageData& image, const LAMBDA& lambda)
{
	size_t pixelIndex = 0;
	for (size_t y = 0; y < image.m_height; ++y)
	{
		SColor* pixel = (SColor*)&image.m_pixels[y * image.m_pitch];
		for (size_t x = 0; x < image.m_width; ++x)
		{
			lambda(*pixel, pixelIndex);
			++pixel;
			++pixelIndex;
		}
	}
}

//======================================================================================
template <typename LAMBDA>
void ImageForEachPixel(const SImageData& image, const LAMBDA& lambda)
{
	size_t pixelIndex = 0;
	for (size_t y = 0; y < image.m_height; ++y)
	{
		SColor* pixel = (SColor*)&image.m_pixels[y * image.m_pitch];
		for (size_t x = 0; x < image.m_width; ++x)
		{
			lambda(*pixel, pixelIndex);
			++pixel;
			++pixelIndex;
		}
	}
}

//======================================================================================
void ImageConvertToLuma(SImageData& image)
{
	ImageForEachPixel(
		image,
		[](SColor& pixel, size_t pixelIndex)
	{
		float luma = float(pixel.R) * 0.3f + float(pixel.G) * 0.59f + float(pixel.B) * 0.11f;
		uint8 lumau8 = uint8(luma + 0.5f);
		pixel.R = lumau8;
		pixel.G = lumau8;
		pixel.B = lumau8;
	}
	);
}

//======================================================================================
void ImageCombine2(const SImageData& imageA, const SImageData& imageB, SImageData& result)
{
	// put the images side by side. A on left, B on right
	ImageInit(result, imageA.m_width + imageB.m_width, max(imageA.m_height, imageB.m_height));
	std::fill(result.m_pixels.begin(), result.m_pixels.end(), 0);

	// image A on left
	for (size_t y = 0; y < imageA.m_height; ++y)
	{
		SColor* destPixel = (SColor*)&result.m_pixels[y * result.m_pitch];
		SColor* srcPixel = (SColor*)&imageA.m_pixels[y * imageA.m_pitch];
		for (size_t x = 0; x < imageA.m_width; ++x)
		{
			destPixel[0] = srcPixel[0];
			++destPixel;
			++srcPixel;
		}
	}

	// image B on right
	for (size_t y = 0; y < imageB.m_height; ++y)
	{
		SColor* destPixel = (SColor*)&result.m_pixels[y * result.m_pitch + imageA.m_width * 3];
		SColor* srcPixel = (SColor*)&imageB.m_pixels[y * imageB.m_pitch];
		for (size_t x = 0; x < imageB.m_width; ++x)
		{
			destPixel[0] = srcPixel[0];
			++destPixel;
			++srcPixel;
		}
	}
}

//======================================================================================
void ImageCombine3(const SImageData& imageA, const SImageData& imageB, const SImageData& imageC, SImageData& result)
{
	// put the images side by side. A on left, B in middle, C on right
	ImageInit(result, imageA.m_width + imageB.m_width + imageC.m_width, max(max(imageA.m_height, imageB.m_height), imageC.m_height));
	std::fill(result.m_pixels.begin(), result.m_pixels.end(), 0);

	// image A on left
	for (size_t y = 0; y < imageA.m_height; ++y)
	{
		SColor* destPixel = (SColor*)&result.m_pixels[y * result.m_pitch];
		SColor* srcPixel = (SColor*)&imageA.m_pixels[y * imageA.m_pitch];
		for (size_t x = 0; x < imageA.m_width; ++x)
		{
			destPixel[0] = srcPixel[0];
			++destPixel;
			++srcPixel;
		}
	}

	// image B in middle
	for (size_t y = 0; y < imageB.m_height; ++y)
	{
		SColor* destPixel = (SColor*)&result.m_pixels[y * result.m_pitch + imageA.m_width * 3];
		SColor* srcPixel = (SColor*)&imageB.m_pixels[y * imageB.m_pitch];
		for (size_t x = 0; x < imageB.m_width; ++x)
		{
			destPixel[0] = srcPixel[0];
			++destPixel;
			++srcPixel;
		}
	}

	// image C on right
	for (size_t y = 0; y < imageC.m_height; ++y)
	{
		SColor* destPixel = (SColor*)&result.m_pixels[y * result.m_pitch + imageA.m_width * 3 + imageC.m_width * 3];
		SColor* srcPixel = (SColor*)&imageC.m_pixels[y * imageC.m_pitch];
		for (size_t x = 0; x < imageC.m_width; ++x)
		{
			destPixel[0] = srcPixel[0];
			++destPixel;
			++srcPixel;
		}
	}
}

//======================================================================================
float GoldenRatioMultiple(size_t multiple)
{
	return float(multiple) * (1.0f + std::sqrtf(5.0f)) / 2.0f;
}

//======================================================================================
void IntegrationTest(const SImageData& dither, const SImageData& groundTruth, size_t frameIndex, const char* label)
{
	// calculate min, max, total and average error
	size_t minError = 0;
	size_t maxError = 0;
	size_t totalError = 0;
	size_t pixelCount = 0;
	for (size_t y = 0; y < dither.m_height; ++y)
	{
		SColor* ditherPixel = (SColor*)&dither.m_pixels[y * dither.m_pitch];
		SColor* truthPixel = (SColor*)&groundTruth.m_pixels[y * groundTruth.m_pitch];
		for (size_t x = 0; x < dither.m_width; ++x)
		{
			size_t error = 0;
			if (ditherPixel->R > truthPixel->R)
				error = ditherPixel->R - truthPixel->R;
			else
				error = truthPixel->R - ditherPixel->R;

			totalError += error;

			if ((x == 0 && y == 0) || error < minError)
				minError = error;

			if ((x == 0 && y == 0) || error > maxError)
				maxError = error;

			++ditherPixel;
			++truthPixel;
			++pixelCount;
		}
	}
	float averageError = float(totalError) / float(pixelCount);

	// calculate standard deviation
	float sumSquaredDiff = 0.0f;
	for (size_t y = 0; y < dither.m_height; ++y)
	{
		SColor* ditherPixel = (SColor*)&dither.m_pixels[y * dither.m_pitch];
		SColor* truthPixel = (SColor*)&groundTruth.m_pixels[y * groundTruth.m_pitch];
		for (size_t x = 0; x < dither.m_width; ++x)
		{
			size_t error = 0;
			if (ditherPixel->R > truthPixel->R)
				error = ditherPixel->R - truthPixel->R;
			else
				error = truthPixel->R - ditherPixel->R;

			float diff = float(error) - averageError;

			sumSquaredDiff += diff*diff;
		}
	}
	float stdDev = std::sqrtf(sumSquaredDiff / float(pixelCount - 1));

	// report results
	fprintf(g_logFile, "%s %zu error\n", label, frameIndex);
	fprintf(g_logFile, "  min error: %zu\n", minError);
	fprintf(g_logFile, "  max error: %zu\n", maxError);
	fprintf(g_logFile, "  avg error: %0.2f\n", averageError);
	fprintf(g_logFile, "  stddev: %0.2f\n", stdDev);
	fprintf(g_logFile, "\n");
}

//======================================================================================
template <size_t NUMFRAMES>
void IntegrationTest2(const SImageData& dither, const SImageData& groundTruth, size_t frameIndex, std::array<size_t, NUMFRAMES>& outMinError, std::array<size_t, NUMFRAMES>& outMaxError, std::array<float, NUMFRAMES>& outAverageError, std::array<float, NUMFRAMES>& outStdDevError)
{
	// calculate min, max, total and average error
	size_t minError = 0;
	size_t maxError = 0;
	size_t totalError = 0;
	size_t pixelCount = 0;
	for (size_t y = 0; y < dither.m_height; ++y)
	{
		SColor* ditherPixel = (SColor*)&dither.m_pixels[y * dither.m_pitch];
		SColor* truthPixel = (SColor*)&groundTruth.m_pixels[y * groundTruth.m_pitch];
		for (size_t x = 0; x < dither.m_width; ++x)
		{
			size_t error = 0;
			if (ditherPixel->R > truthPixel->R)
				error = ditherPixel->R - truthPixel->R;
			else
				error = truthPixel->R - ditherPixel->R;

			totalError += error;

			if ((x == 0 && y == 0) || error < minError)
				minError = error;

			if ((x == 0 && y == 0) || error > maxError)
				maxError = error;

			++ditherPixel;
			++truthPixel;
			++pixelCount;
		}
	}
	float averageError = float(totalError) / float(pixelCount);

	// calculate standard deviation
	float sumSquaredDiff = 0.0f;
	for (size_t y = 0; y < dither.m_height; ++y)
	{
		SColor* ditherPixel = (SColor*)&dither.m_pixels[y * dither.m_pitch];
		SColor* truthPixel = (SColor*)&groundTruth.m_pixels[y * groundTruth.m_pitch];
		for (size_t x = 0; x < dither.m_width; ++x)
		{
			size_t error = 0;
			if (ditherPixel->R > truthPixel->R)
				error = ditherPixel->R - truthPixel->R;
			else
				error = truthPixel->R - ditherPixel->R;

			float diff = float(error) - averageError;

			sumSquaredDiff += diff*diff;
		}
	}
	float stdDev = std::sqrtf(sumSquaredDiff / float(pixelCount - 1));

	outMinError[frameIndex] = minError;
	outMaxError[frameIndex] = maxError;
	outAverageError[frameIndex] = averageError;
	outStdDevError[frameIndex] = stdDev;
}

//======================================================================================
template <size_t NUMFRAMES>
void WriteIntegrationTest2(std::array<size_t, NUMFRAMES>& minError, std::array<size_t, NUMFRAMES>& maxError, std::array<float, NUMFRAMES>& averageError, std::array<float, NUMFRAMES>& stdDevError, const char* label)
{
	fprintf(g_logFile, "%s error\n", label);
	fprintf(g_logFile, "  min error: ");
	for (size_t v : minError)
		fprintf(g_logFile, "%zu, ", v);

	fprintf(g_logFile, "\n  max error: ");
	for (size_t v : maxError)
		fprintf(g_logFile, "%zu, ", v);

	fprintf(g_logFile, "\n  avg error: ");
	for (float v : averageError)
		fprintf(g_logFile, "%0.2f, ", v);

	fprintf(g_logFile, "\n  stddev: ");
	for (float v : stdDevError)
		fprintf(g_logFile, "%0.2f, ", v);

	fprintf(g_logFile, "\n\n");
}

//======================================================================================
void HistogramTest(const SImageData& noise, size_t frameIndex, const char* label)
{
	std::array<size_t, 256> counts;
	std::fill(counts.begin(), counts.end(), 0);

	ImageForEachPixel(
		noise,
		[&](const SColor& pixel, size_t pixelIndex)
	{
		counts[pixel.R]++;
	}
	);

	// calculate min, max, total and average
	size_t minCount = 0;
	size_t maxCount = 0;
	size_t totalCount = 0;
	for (size_t i = 0; i < 256; ++i)
	{
		if (i == 0 || counts[i] < minCount)
			minCount = counts[i];

		if (i == 0 || counts[i] > maxCount)
			maxCount = counts[i];

		totalCount += counts[i];
	}
	float averageCount = float(totalCount) / float(256.0f);

	// calculate standard deviation
	float sumSquaredDiff = 0.0f;
	for (size_t i = 0; i < 256; ++i)
	{
		float diff = float(counts[i]) - averageCount;
		sumSquaredDiff += diff*diff;
	}
	float stdDev = std::sqrtf(sumSquaredDiff / 255.0f);

	// report results
	fprintf(g_logFile, "%s %zu histogram\n", label, frameIndex);
	fprintf(g_logFile, "  min count: %zu\n", minCount);
	fprintf(g_logFile, "  max count: %zu\n", maxCount);
	fprintf(g_logFile, "  avg count: %0.2f\n", averageCount);
	fprintf(g_logFile, "  stddev: %0.2f\n", stdDev);
	fprintf(g_logFile, "  counts: ");
	for (size_t i = 0; i < 256; ++i)
	{
		if (i > 0)
			fprintf(g_logFile, ", ");
		fprintf(g_logFile, "%zu", counts[i]);
	}

	fprintf(g_logFile, "\n\n");
}

//======================================================================================
void GenerateWhiteNoise(SImageData& image, size_t width, size_t height)
{
	ImageInit(image, width, height);

	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<unsigned int> dist(0, 255);

	ImageForEachPixel(
		image,
		[&](SColor& pixel, size_t pixelIndex)
	{
		uint8 value = dist(rng);
		pixel.R = value;
		pixel.G = value;
		pixel.B = value;
	}
	);
}

//======================================================================================
void GenerateInterleavedGradientNoise(SImageData& image, size_t width, size_t height, float offsetX, float offsetY)
{
	ImageInit(image, width, height);

	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<unsigned int> dist(0, 255);

	for (size_t y = 0; y < height; ++y)
	{
		SColor* pixel = (SColor*)&image.m_pixels[y * image.m_pitch];
		for (size_t x = 0; x < width; ++x)
		{
			float valueFloat = std::fmodf(52.9829189f * std::fmod(0.06711056f*float(x + offsetX) + 0.00583715f*float(y + offsetY), 1.0f), 1.0f);
			size_t valueBig = size_t(valueFloat * 256.0f);
			uint8 value = uint8(valueBig % 256);
			pixel->R = value;
			pixel->G = value;
			pixel->B = value;
			++pixel;
		}
	}
}

//======================================================================================
template <size_t NUM_SAMPLES>
void GenerateVanDerCoruptSequence(std::array<float, NUM_SAMPLES>& samples, size_t base)
{
	for (size_t i = 0; i < NUM_SAMPLES; ++i)
	{
		samples[i] = 0.0f;
		float denominator = float(base);
		size_t n = i;
		while (n > 0)
		{
			size_t multiplier = n % base;
			samples[i] += float(multiplier) / denominator;
			n = n / base;
			denominator *= base;
		}
	}
}

//======================================================================================
void DitherWithTexture(const SImageData& ditherImage, const SImageData& noiseImage, SImageData& result)
{
	// init the result image
	ImageInit(result, ditherImage.m_width, ditherImage.m_height);

	// make the result image
	for (size_t y = 0; y < ditherImage.m_height; ++y)
	{
		SColor* srcDitherPixel = (SColor*)&ditherImage.m_pixels[y * ditherImage.m_pitch];
		SColor* destDitherPixel = (SColor*)&result.m_pixels[y * result.m_pitch];

		for (size_t x = 0; x < ditherImage.m_width; ++x)
		{
			// tile the noise in case it isn't the same size as the image we are dithering
			size_t noiseX = x % noiseImage.m_width;
			size_t noiseY = y % noiseImage.m_height;
			SColor* noisePixel = (SColor*)&noiseImage.m_pixels[noiseY * noiseImage.m_pitch + noiseX * 3];

			uint8 value = 0;
			if (noisePixel->R < srcDitherPixel->R)
				value = 255;

			destDitherPixel->R = value;
			destDitherPixel->G = value;
			destDitherPixel->B = value;

			++srcDitherPixel;
			++destDitherPixel;
		}
	}
}
//======================================================================================
void DitherInterleavedGradientNoise(const SImageData& ditherImage)
{
	printf("\n%s\n", __FUNCTION__);

	// make noise
	SImageData noise;
	GenerateInterleavedGradientNoise(noise, ditherImage.m_width, ditherImage.m_height, 0.0f, 0.0f);

	// dither the image
	SImageData dither;
	DitherWithTexture(ditherImage, noise, dither);


	SImageDataComplex noiseDFT;
	SImageData noiseDFTMag;

	//if (c_doDFT)
	//{
	ImageDFT(dither, noiseDFT);
	GetMagnitudeData(noiseDFT, noiseDFTMag);
	//}
	//else
	//{
	//	ImageInit(noiseDFTMag, ivNoise.m_width, ivNoise.m_height);
	//	std::fill(noiseDFTMag.m_pixels.begin(), noiseDFTMag.m_pixels.end(), 0);
	//}

	// Histogram test the noise
	//	HistogramTest(ivNoise, i, __FUNCTION__);


	//DitherWithTexture(ditherImage, ivNoise, dither);

	// save the results
	SImageData combined;
	ImageCombine3(noise, noiseDFTMag, dither, combined);
	ImageSave(combined, "out/still_ignoise.bmp");
}

//======================================================================================
void DitherBlueNoise(const SImageData& ditherImage, const SImageData& blueNoise)
{
	printf("\n%s\n", __FUNCTION__);

	// dither the image
	SImageData dither;
	DitherWithTexture(ditherImage, blueNoise, dither);


	SImageDataComplex noiseDFT;
	SImageData noiseDFTMag;

	//if (c_doDFT)
	//{
	ImageDFT(dither, noiseDFT);
	GetMagnitudeData(noiseDFT, noiseDFTMag);
	//}
	//else
	//{
	//	ImageInit(noiseDFTMag, ivNoise.m_width, ivNoise.m_height);
	//	std::fill(noiseDFTMag.m_pixels.begin(), noiseDFTMag.m_pixels.end(), 0);
	//}

	// Histogram test the noise
	//	HistogramTest(ivNoise, i, __FUNCTION__);


	//DitherWithTexture(ditherImage, ivNoise, dither);

	// save the results
	SImageData combined;
	ImageCombine3(blueNoise, noiseDFTMag, dither, combined);

	ImageSave(combined, "out/still_bluenoise.bmp");
}

void ImageResize(SImageData &img_orig, SImageData &img_resized, int target_width, int target_height)
{
	ImageInit(img_resized, target_width, target_height);



	float X, Y, x, y, a, b;

	for (int i = 0; i < img_resized.m_width; i++)
		for (int j = 0; j < img_resized.m_height; j++)
		{
			SColor & targetC = img_resized.GetPixel(i, j);

			x = img_orig.m_width * (1.0 * i / img_resized.m_width);	//Scaling in x direction
			y = img_orig.m_height * (1.0 * j / img_resized.m_height);	//Scaling in x direction

			X = floor(x);
			Y = floor(y);

			SColor & c1 = img_orig.GetPixel(X, Y);// *(1.0 - a)*(1.0 - b);
			SColor & c2 = img_orig.GetPixel(X, Y + 1);// *(1.0 - a)*(1.0 - b);
			SColor & c3 = img_orig.GetPixel(X + 1, Y);// *(1.0 - a)*(1.0 - b);
			SColor & c4 = img_orig.GetPixel(X + 1, Y + 1);// *(1.0 - a)*(1.0 - b);


			targetC.R = c1.R * (1.0 - a)*(1.0 - b) + c2.R * (a)*(1.0 - b) + c3.R * (1.0 - a)*(b)+c4.R * a *b;
			targetC.G = c1.G * (1.0 - a)*(1.0 - b) + c2.G * (a)*(1.0 - b) + c3.G * (1.0 - a)*(b)+c4.G * a *b;
			targetC.B = c1.B * (1.0 - a)*(1.0 - b) + c2.B * (a)*(1.0 - b) + c3.B * (1.0 - a)*(b)+c4.B * a *b;
			//BILINEAR INTERPOLATION FORMULATION
		}
}
// Image Resize without BILINEAR INTERPOLATION 
void ImageResizeNI(SImageData &img_orig, SImageData &img_resized, int target_width, int target_height)
{
	ImageInit(img_resized, target_width, target_height);



	float X, Y, x, y, a, b;

	for (int i = 0; i < img_resized.m_width; i++)
		for (int j = 0; j < img_resized.m_height; j++)
		{
			SColor & targetC = img_resized.GetPixel(i, j);

			x = img_orig.m_width * (1.0 * i / img_resized.m_width);	//Scaling in x direction
			y = img_orig.m_height * (1.0 * j / img_resized.m_height);	//Scaling in x direction

			X = floor(x);
			Y = floor(y);

			a = y - Y;
			b = x - X;

			SColor & c1 = img_orig.GetPixel(X, Y);// *(1.0 - a)*(1.0 - b);

			targetC.R = c1.R;// *(1.0 - a)*(1.0 - b) + c2.R * (a)*(1.0 - b) + c3.R * (1.0 - a)*(b)+c4.R * a *b;
			targetC.G = c1.G;// *(1.0 - a)*(1.0 - b) + c2.G * (a)*(1.0 - b) + c3.G * (1.0 - a)*(b)+c4.G * a *b;
			targetC.B = c1.B;// *(1.0 - a)*(1.0 - b) + c2.B * (a)*(1.0 - b) + c3.B * (1.0 - a)*(b)+c4.B * a *b;
			//BILINEAR INTERPOLATION FORMULATION
		}
}

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

struct simpleFPS
{
	double x;
	double y;
};
std::vector<simpleFPS> my_pre_computed_samples;

bool LoadFPSSeqs(const char * name, int target_size)
{
	my_pre_computed_samples.clear();
	std::ifstream myfile;
	myfile.open(name, ios::in);
	if (myfile.is_open())
	{
		long size = 0;
		myfile >> size;
		if (target_size < size)
			size = target_size;
		my_pre_computed_samples.resize(size);
		for (long i = 0; i < size; ++i)
		{
			myfile >> my_pre_computed_samples[i].x >> my_pre_computed_samples[i].y;
			//my_pre_computed_samples[i].seqs = i; // TO FIX: ²»±ØÒª£¿
		}
		std::cout << "load " << name << "  success!" << std::endl;
		std::cout << "total size: " << my_pre_computed_samples.size() << "!" << std::endl;
		long i = my_pre_computed_samples.size() - 1;
		std::cout << "Last value: " << my_pre_computed_samples[i].x << " " << my_pre_computed_samples[i].y << std::endl;
		myfile.close();

		return true;
	}
	else
	{
		std::cout << "load " << name << "  failed!" << std::endl;
		return false;
	}
}
//======================================================================================


#include "_fps.h"

extern double fps[524288][2];

void GenerateFPSNoise(SImageData& image, size_t width, size_t height)
{
	ImageInit(image, width, height);

	float overSc = 0.5;

	for (int i = 0; i < width * height * overSc; ++i)
	{
		double x = fps[i][0];
		double y = fps[i][1];

		size_t k = x * width;
		size_t t = y * height;

		if (k < 0 || k >= width || t <0 || t> height)
		{
			cout << "error int fps data\n";
			continue;
		}
		//int32(out(i, j) * 255 / (msize * msize));
		uint8 value = (i + 1) * 255 / (width * height * overSc);
		SColor& pixel = image.GetPixel(k, t);
		pixel.R = value;
		pixel.G = value;
		pixel.B = value;

	}
}

//======================================================================================
void DitherIVNoise(const SImageData& ditherImage, bool doDFT = true, int num = 0)
{
	printf("\n%s\n", __FUNCTION__);

	SImageData dither;
	SImageData dither2;
	SImageData noise;
	SImageData noise2;
	ImageInit(noise, ditherImage.m_width, ditherImage.m_height);
	ImageInit(dither, ditherImage.m_width, ditherImage.m_height);
	ImageInit(dither2, ditherImage.m_width * 2, ditherImage.m_height * 2);
	ImageInit(noise2, ditherImage.m_width * 2, ditherImage.m_height * 2);
	//SImageData dither2;
	size_t width = ditherImage.m_width;
	size_t height = ditherImage.m_height;
	// make noise
	float overSc = 1.10;
	float _f = 1.0 / (width * height * overSc);

	int sampleIdx = 1;
	overSc = 1.1;
	sampleIdx = 1;
	for (int i = 0; i < (width * height * overSc); ++i)
	{
		double x = fps[i][0];
		double y = fps[i][1];

		size_t k = x * width;
		size_t t = y * height;


		size_t k2 = x * width * 2;
		size_t t2 = y * height * 2;

		if (k < 0 || k >= width || t <0 || t> height)
		{
			cout << "error int fps data\n";
			continue;
		}

		uint8 value = (sampleIdx) * 255 / (width * height * 1);		

		SColor& noisePixel = noise2.GetPixel(k, t);
		const SColor& pixel = ditherImage.GetPixel(k, t);
		SColor& destDitherPixel = dither.GetPixel(k, t);
		SColor& destDither2Pixel = dither2.GetPixel(k2, t2);

		if (noisePixel.R == 255)
		{
			noisePixel.R = int(value);// / 16) * 16;
			noisePixel.G = int(value);// / 16) * 16;
			noisePixel.B = int(value);// / 16) * 16;
			sampleIdx++;

			uint8 R = 255 - pixel.R;

			uint8 Rvalue = 255;
			if (R > value)
				Rvalue = 0;

			destDither2Pixel.R = Rvalue;
			destDither2Pixel.G = Rvalue;
			destDither2Pixel.B = Rvalue;


			destDitherPixel.R = Rvalue;
			destDitherPixel.G = Rvalue;
			destDitherPixel.B = Rvalue;
		}


	}

	char buffer[256];
	sprintf(buffer, "out/IVDithered%d.bmp", num);

	ImageSave(dither, buffer);
	sprintf(buffer, "out/IVDithered2_%d.bmp", num);
	ImageSave(dither2, buffer);

	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<unsigned int> dist(0, 255);
	uint8 value = dist(rng);

	if (doDFT)
	{
		SImageDataComplex noiseDFT;
		SImageData noiseDFTMag;
		ImageDFT(dither2, noiseDFT);

		GetMagnitudeData(noiseDFT, noiseDFTMag);
		

		sprintf(buffer, "out/IVDithered2_DFT%d.bmp", num);
		ImageSave(noiseDFTMag, buffer);
	}
}


//======================================================================================
void DitherIVNoise(const SImageData& ditherImage, const char * file_size_type_Name, bool doDFT = true)
{
	printf("\n%s\n", __FUNCTION__);

	std::vector<int64> pixelDifferences;
	pixelDifferences.resize(ditherImage.m_width * ditherImage.m_height);


	SImageData dither;
	SImageData dither2;
	SImageData noise;
	SImageData noise2;
	ImageInit(noise, ditherImage.m_width, ditherImage.m_height);
	ImageInit(dither, ditherImage.m_width, ditherImage.m_height);
	ImageInit(dither2, ditherImage.m_width * 2, ditherImage.m_height * 2);
	ImageInit(noise2, ditherImage.m_width * 2, ditherImage.m_height * 2);
	//SImageData dither2;
	size_t width = ditherImage.m_width;
	size_t height = ditherImage.m_height;
	// make noise
	float overSc = 1.50;
	float _f = 1.0 / (width * height * overSc);

	int sampleIdx = 1;
	overSc = 1.1;
	sampleIdx = 1;
	for (int i = 0; i < (width * height * overSc); ++i)
	{
		double x = fps[i][0];
		double y = fps[i][1];

		size_t k = x * width;
		size_t t = y * height;


		size_t k2 = x * width * 2;
		size_t t2 = y * height * 2;

		if (k < 0 || k >= width || t <0 || t> height)
		{
			cout << "error int fps data\n";
			continue;
		}

		uint8 value = (sampleIdx) * 255 / (width * height );

		SColor& noisePixel = noise.GetPixel(k, t);
		const SColor& pixel = ditherImage.GetPixel(k, t);
		SColor& destDitherPixel = dither.GetPixel(k, t);
		SColor& destDither2Pixel = dither2.GetPixel(k2, t2);

		if (noisePixel.R == 255)
		{
			noisePixel.R = int(value);// / 16) * 16;
			noisePixel.G = int(value);// / 16) * 16;
			noisePixel.B = int(value);// / 16) * 16;
			sampleIdx++;

			uint8 R = 255 - pixel.R;

			uint8 Rvalue = 255;
			if (R > value)
				Rvalue = 0;

			

			destDither2Pixel.R = Rvalue;
			destDither2Pixel.G = Rvalue;
			destDither2Pixel.B = Rvalue;


			destDitherPixel.R = Rvalue;
			destDitherPixel.G = Rvalue;
			destDitherPixel.B = Rvalue;

			pixelDifferences[t * ditherImage.m_width + k] = int64(destDitherPixel.R) - int64(pixel.R);

		}


	}


	char outFileName[256];
	strcpy(outFileName, file_size_type_Name);
	strcat(outFileName, "_noise.bmp");
	ImageSave(noise, outFileName);
	
	strcpy(outFileName, file_size_type_Name);
	strcat(outFileName, "_stipple.bmp");
	ImageSave(dither, outFileName);

	strcpy(outFileName, file_size_type_Name);
	strcat(outFileName, "_stipple2.bmp");
	ImageSave(dither2, outFileName);



	// calculate some metrics
	int64 totalDiff = 0;
	int64 totalAbsDiff = 0;
	for (int64 v : pixelDifferences)
	{
		totalDiff += v;
		totalAbsDiff += std::abs(v);
	}
	float averageDiff = float(totalDiff) / float(pixelDifferences.size());
	float stdDev = 0.0f;
	for (int64 v : pixelDifferences)
	{
		stdDev += (float(v) - averageDiff) * (float(v) - averageDiff);
	}
	stdDev = std::sqrt(stdDev / float(pixelDifferences.size()));

	fprintf(g_logFile, "%s\nTotal Diff: %zi\nTotal Abs Diff: %zi\nAverage Diff:%0.2f\nStdDev. Diff: %0.2f\n\n", file_size_type_Name, totalDiff, totalAbsDiff, averageDiff, stdDev);



	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<unsigned int> dist(0, 255);
	uint8 value = dist(rng);

	if (doDFT)
	{
		SImageDataComplex noiseDFT;
		SImageData noiseDFTMag;
		ImageDFT(dither2, noiseDFT);

		GetMagnitudeData(noiseDFT, noiseDFTMag);
		SImageData noiseDFTMagResized;
		ImageResize(noiseDFTMag, noiseDFTMagResized, noiseDFTMag.m_width / 2, noiseDFTMag.m_height / 2);

		char buffer[256];
		strcpy(buffer, file_size_type_Name);
		strcat(buffer, "_dftMag.bmp");
		ImageSave(noiseDFTMagResized, buffer);
	}
}

//======================================================================================
void DitherWithNoiseTexture(const SImageData& ditherImage, const char * outputStr, const SImageData& noiseImage, SImageData& result)
{
	// init the result image
	ImageInit(result, ditherImage.m_width, ditherImage.m_height);

	std::vector<int64> pixelDifferences;
	pixelDifferences.resize(ditherImage.m_width * ditherImage.m_height);

	// make the result image
	for (size_t y = 0; y < ditherImage.m_height; ++y)
	{
		SColor* srcDitherPixel = (SColor*)&ditherImage.m_pixels[y * ditherImage.m_pitch];
		SColor* destDitherPixel = (SColor*)&result.m_pixels[y * result.m_pitch];

		for (size_t x = 0; x < ditherImage.m_width; ++x)
		{
			// tile the noise in case it isn't the same size as the image we are dithering
			size_t noiseX = x % noiseImage.m_width;
			size_t noiseY = y % noiseImage.m_height;
			SColor* noisePixel = (SColor*)&noiseImage.m_pixels[noiseY * noiseImage.m_pitch + noiseX * 3];

			
			if (noisePixel->R >= srcDitherPixel->R)				
				destDitherPixel->Set(0, 0, 0);
			else
				destDitherPixel->Set(255, 255, 255);

			pixelDifferences[y * ditherImage.m_width + x] = int64(destDitherPixel->R) - int64(srcDitherPixel->R);

			++srcDitherPixel;
			++destDitherPixel;			
		}
	}


	// calculate some metrics
	int64 totalDiff = 0;
	int64 totalAbsDiff = 0;
	for (int64 v : pixelDifferences)
	{
		totalDiff += v;
		totalAbsDiff += std::abs(v);
	}
	float averageDiff = float(totalDiff) / float(pixelDifferences.size());
	float stdDev = 0.0f;
	for (int64 v : pixelDifferences)
	{
		stdDev += (float(v) - averageDiff) * (float(v) - averageDiff);
	}
	stdDev = std::sqrt(stdDev / float(pixelDifferences.size()));

	fprintf(g_logFile, "%s\nTotal Diff: %zi\nTotal Abs Diff: %zi\nAverage Diff:%0.2f\nStdDev. Diff: %0.2f\n\n", outputStr, totalDiff, totalAbsDiff, averageDiff, stdDev);

	// save the image
	char outFileName[256];
	strcpy(outFileName, outputStr);
	strcat(outFileName, "_stipple.bmp");
	ImageSave(result, outFileName);
}


//======================================================================================
void DitherInputNoise(const SImageData& ditherImage, const SImageData& inputNoise,const char * file_size_type_Name,bool doDFT = true)
{
	printf("\n%s\n", __FUNCTION__);
	// dither the image

	SImageData dither;
	DitherWithNoiseTexture(ditherImage, file_size_type_Name,inputNoise, dither);
	
	if (doDFT)
	{
		SImageDataComplex noiseDFT;
		SImageData noiseDFTMag;
		ImageDFT(dither, noiseDFT);
		GetMagnitudeData(noiseDFT, noiseDFTMag);
		char outFileName[256];
		strcpy(outFileName, file_size_type_Name);
		strcat(outFileName, "_dftMag.bmp");

		ImageSave(noiseDFTMag, outFileName);
	}
}




void DitherAllWithDegrateResolution(SImageData &ditherImage ,const char * fileName)
{
	//
	//cout << "Step 1: Resized Image\n";
	//cout << "Step 1: Resized Image\n";
	char ditherfileName[256];
	char num[10];
	bool doDFT = false;

	for (int i = 7; i <= 8; ++i)
	{
		if (i >= 7 && i <= 8)
			doDFT = false;
		else 
			doDFT = false;

		char NoiseBuffer[256];
		int w = pow(2, i);
		SImageData img_resized;
	
		strcpy(ditherfileName, fileName);
		sprintf(num, "_%d", w);
		strcat(ditherfileName, num);
		
		char outFileName[256];
		strcpy(outFileName, ditherfileName);
		strcat(outFileName, ".bmp");

		ImageResize(ditherImage, img_resized,w, w);		
		ImageSave(img_resized, outFileName);
		
		//char file_size_type_Name[256];
		///*
		//	Dither with bayer matrix.
		//*/
		//SImageData BayerNoise;
		//sprintf(NoiseBuffer, "src/Bayer_%d.bmp", w);
		//if (!ImageLoad(NoiseBuffer, BayerNoise))
		//{
		//	printf("Could not load %s", NoiseBuffer);
		//	return;
		//}
		//
		//strcpy(file_size_type_Name, ditherfileName);
		//strcat(file_size_type_Name, "Bayer");
		//DitherInputNoise(img_resized, BayerNoise, file_size_type_Name, doDFT);

		///*
		//	Dither with VNC mask.
		//*/
		//SImageData VNCNoise;		
		//sprintf(NoiseBuffer, "src/VNC%d.bmp",w);
		//if (!ImageLoad(NoiseBuffer, VNCNoise))
		//{
		//	printf("Could not load %s", NoiseBuffer);
		//	return ;
		//}
		//strcpy(file_size_type_Name, ditherfileName);
		//strcat(file_size_type_Name, "VNC");		
		//DitherInputNoise(img_resized, VNCNoise, file_size_type_Name, doDFT);
		//

		///*
		//	Dither with blue 16 mask.
		//*/
		//SImageData blueNoise;
		//sprintf(NoiseBuffer, "src/BN%d.bmp", w);
		//if (!ImageLoad(NoiseBuffer, blueNoise))
		//{
		//	printf("Could not load %s", NoiseBuffer);
		//	return;
		//}

		//strcpy(file_size_type_Name, ditherfileName);
		//strcat(file_size_type_Name, "blue16");
		//DitherInputNoise(img_resized, blueNoise, file_size_type_Name, doDFT);


		///*
		//Dither with blue high pass mask.
		//*/
		//SImageData blueNoiseHP;
		//sprintf(NoiseBuffer, "src/BNHP%d.bmp", w);
		//if (!ImageLoad(NoiseBuffer, blueNoiseHP))
		//{
		//	printf("Could not load %s", NoiseBuffer);
		//	return;
		//}

		//strcpy(file_size_type_Name, ditherfileName);
		//strcat(file_size_type_Name, "blueHighPass");
		//DitherInputNoise(img_resized, blueNoiseHP, file_size_type_Name, doDFT);
		//	
		//strcpy(file_size_type_Name, ditherfileName);
		//strcat(file_size_type_Name, "_RIIV_");

		//DitherIVNoise(img_resized, file_size_type_Name, doDFT);

		
	}
}

#include<direct.h>
//======================================================================================
// Load 
void _Main(string file) 
{
	std::string filename;
	filename = file;

	mkdir(filename.c_str());


	string logFile = filename + "_log.txt";
	g_logFile = fopen(logFile.c_str(), "w+t");

	string realFileName = filename + ".bmp";
	printf("Try to load %s\n", filename.c_str());

	string filePath = filename + "/" + filename;

	SImageData ditherImage;
	if (!ImageLoad(realFileName.c_str(), ditherImage))
	{
		printf("Could not load %s", realFileName.c_str());
			return;	

	}

	ImageConvertToLuma(ditherImage);
	DitherAllWithDegrateResolution(ditherImage, filePath.c_str());
	fclose(g_logFile);

}

string fileList[] = 
{
	"gray240",
	"gradient",
	"Face",
	"BoneChina",
	"BabyHand",
	"shoe",
	"plant"
};


//======================================================================================
void ResolutionFreeIVNoiseTestWithLinearGradient(const char * prefix, int width, int height, bool doDFT = true)
{
	printf("\n%s\n", __FUNCTION__);

	std::vector<int64> pixelDifferences;
	pixelDifferences.resize(width * height);
		

	SImageData ditherGray;
	SImageData dithergradient;
	SImageData noise, noise2;

	
	ImageInit(ditherGray, width, height);
	ImageInit(dithergradient, width, height);
	ImageInit(noise, width, height);
	ImageInit(noise2, width * 2, height * 2);

	// make noise
	float overSc = 1.50;
	float _f = 1.0 / (width * height * overSc);

	int sampleIdx = 1;
	overSc = 1.1;
	sampleIdx = 1;
	for (int i = 0; i < (width * height * overSc); ++i)
	{
		double x = fps[i][0];
		double y = fps[i][1];

		size_t k = x * width;
		size_t t = y * height;


		size_t k2 = x * width * 2;
		size_t t2 = y * height * 2;

		if (k < 0 || k >= width || t <0 || t> height)
		{
			cout << "error int fps data\n";
			continue;
		}

		uint8 value = (sampleIdx) * 255 / (width * height);

		
		const SColor pixel(240,240,240);

		SColor& noisePixel = noise.GetPixel(k, t);
		//SColor& noisePixel = noise.GetPixel(k2, t2);
		SColor& dstDitherGrayPixel = ditherGray.GetPixel(k, t);
		SColor& dstDitherGradientPixel = dithergradient.GetPixel(k, t);

		float gradientPixel = 255.0f * float(k) / float(width);

		

		if (noisePixel.R == 255)
		{
			noisePixel.R = int(value);// / 16) * 16;
			noisePixel.G = int(value);// / 16) * 16;
			noisePixel.B = int(value);// / 16) * 16;
			sampleIdx++;

			
			if (15 > value) 
			{
				dstDitherGrayPixel.R = 0;
				dstDitherGrayPixel.G = 0;
				dstDitherGrayPixel.B = 0;
			}

			if (gradientPixel > value) 
			{
				dstDitherGradientPixel.R = 0;
				dstDitherGradientPixel.G = 0;
				dstDitherGradientPixel.B = 0;
			}

		//	pixelDifferences[t * ditherImage.m_width + k] = int64(destDitherPixel.R) - int64(pixel.R);

		}


	}


	char outFileName[256];
	strcpy(outFileName, prefix);
	strcat(outFileName, "_noise.bmp");
	ImageSave(noise, outFileName);

	strcpy(outFileName, prefix);
	strcat(outFileName, "_constant.bmp");
	ImageSave(ditherGray, outFileName);

	strcpy(outFileName, prefix);
	strcat(outFileName, "_gradient.bmp");
	ImageSave(dithergradient, outFileName);



	//// calculate some metrics
	//int64 totalDiff = 0;
	//int64 totalAbsDiff = 0;
	//for (int64 v : pixelDifferences)
	//{
	//	totalDiff += v;
	//	totalAbsDiff += std::abs(v);
	//}
	//float averageDiff = float(totalDiff) / float(pixelDifferences.size());
	//float stdDev = 0.0f;
	//for (int64 v : pixelDifferences)
	//{
	//	stdDev += (float(v) - averageDiff) * (float(v) - averageDiff);
	//}
	//stdDev = std::sqrt(stdDev / float(pixelDifferences.size()));

	//fprintf(g_logFile, "%s\nTotal Diff: %zi\nTotal Abs Diff: %zi\nAverage Diff:%0.2f\nStdDev. Diff: %0.2f\n\n", file_size_type_Name, totalDiff, totalAbsDiff, averageDiff, stdDev);



	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<unsigned int> dist(0, 255);
	uint8 value = dist(rng);

	if (doDFT)
	{
		SImageDataComplex noiseDFT;
		SImageData noiseDFTMag;
		ImageDFT(ditherGray, noiseDFT);

		GetMagnitudeData(noiseDFT, noiseDFTMag);
	//	SImageData noiseDFTMagResized;
	//	ImageResize(noiseDFTMag, noiseDFTMagResized, noiseDFTMag.m_width / 2, noiseDFTMag.m_height / 2);

		char buffer[256];
		strcpy(buffer, prefix);
		strcat(buffer, "_dftMag.bmp");
		ImageSave(noiseDFTMag, buffer);
	}
}

// Save IV Noise Mask and related Dithered Result to a directory.
void fun2(string file)
{
	std::string filename;
	filename = file;

	mkdir(filename.c_str());
	string filePath = filename + "/" + filename;
	char Num[256];
	
	for (int i = 4; i < 104; ++i)
	{	
		sprintf(Num, "%dx%d.bmp", i * 5,i);
		ResolutionFreeIVNoiseTestWithLinearGradient((filePath + Num).c_str(), i * 5, i, false);
	}

}

// dithering input images from the filelist with different masks (different resolutions)
int fun1(int argc, char** argv)
{

	if (argc > 1)
	{
		_Main(argv[1]);
	}
	else
	{
		for (auto i : fileList)
			_Main(i);
	}
	return 1;
}

int main(int argc, char** argv)
{
	return fun1(argc, argv);

	//fun2("IVNoise");

}

/*

*/