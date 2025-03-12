#include "helperFunctions.h"
#include <iostream>
#include <bitset>

float* boundingBox = new float[7];
int* numNeighbors = nullptr;

void boundingBoxConstruction(Particle* particles, int numParticles, float h) // TODO reduce calculations? if grid does not change
{
    // xmin, xmax, ymin, ymax, cellsx, cellsy
    boundingBox[0] = particles[0].position.x; // min
    boundingBox[1] = particles[0].position.x; // max
    boundingBox[2] = particles[0].position.y; // min
    boundingBox[3] = particles[0].position.y; // max

    // compute bounding box
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].position.x < boundingBox[0])
        {
            boundingBox[0] = particles[i].position.x;
        }
        else if (particles[i].position.x > boundingBox[1])
        {
            boundingBox[1] = particles[i].position.x;
        }
        if (particles[i].position.y < boundingBox[2])
        {
            boundingBox[2] = particles[i].position.y;
        }
        else if (particles[i].position.y > boundingBox[3])
        {
            boundingBox[3] = particles[i].position.y;
        }
    }

    boundingBox[0] -= 0.4 * h; // make grid slightly bigger for edge cases
    boundingBox[1] += 0.4 * h;
    boundingBox[2] -= 0.4 * h;
    boundingBox[3] += 0.4 * h;

    boundingBox[4] = ceil((boundingBox[1] - boundingBox[0]) / (2.f * h));
    boundingBox[5] = ceil((boundingBox[3] - boundingBox[2]) / (2.f * h));
    boundingBox[6] = boundingBox[4] * boundingBox[5];
}

unsigned int hashFunction(int cellIndexX, int cellIndexY, unsigned int sizeHashTable)
{
    return ((cellIndexX * 73856093u) ^ (cellIndexY * 19349663u)) % sizeHashTable; // overflow creates negative number, either add sizeHashTable, use unsigned int or other larger datatypes
}

// Helper function to spread the bits of a 32-bit integer into 64 bits // https://lemire.me/blog/2018/01/08/how-fast-can-you-bit-interleave-32-bit-integers/
uint64_t spreadBits(uint32_t x) {
    uint64_t result = x;
    result = (result | (result << 16)) & 0x0000FFFF0000FFFF;
    result = (result | (result << 8)) & 0x00FF00FF00FF00FF;
    result = (result | (result << 4)) & 0x0F0F0F0F0F0F0F0F;
    result = (result | (result << 2)) & 0x3333333333333333;
    result = (result | (result << 1)) & 0x5555555555555555;
    return result;
}

// Function to interleave bits of two 32-bit integers x and y
uint64_t interleaveBits(uint32_t x, uint32_t y) {
    return (spreadBits(x) | (spreadBits(y) << 1));
}

bool comp(int a, int b) {
    return a < b;
}

vector<unsigned char> compress(vector<int> particles)
{
    sort(particles.begin(), particles.end(), comp);

    vector<unsigned char> byteStream;
    byteStream.reserve(4);

    unsigned char firstElement[4];
    memcpy(firstElement, &particles[0], sizeof(int));
    for (int j = 3; j >= 0; --j) {
        byteStream.push_back(firstElement[j]);
    }

    unsigned char packedByte = 0;
    int shiftAmount = 6;

    vector<unsigned char> data;

    for (size_t i = 1; i < particles.size(); ++i) {
        unsigned char value = 0;
        int delta = particles[i] - particles[i - 1] - 1;

        if (delta == 0) {
            value = 0b00;  // 00
        }
        else if (delta == 1) {
            value = 0b01;  // 01
        }
        else if (delta >= 2 && delta < 256) {
            value = 0b10;  // 10
            data.push_back(static_cast<unsigned char>(delta));
        }
        else if (delta >= 256) {
            value = 0b11;  // 11
            unsigned char deltaBytes[4];
            memcpy(deltaBytes, &delta, sizeof(int));
            for (int j = 3; j >= 0; --j) {
                data.push_back(deltaBytes[j]);
            }
        }

        packedByte |= (value << shiftAmount);
        shiftAmount -= 2;

        if (shiftAmount < 0) {
            byteStream.push_back(packedByte);
            packedByte = 0;
            shiftAmount = 6;
        }
    }

    if (shiftAmount < 6) {
        byteStream.push_back(packedByte);
    }

    byteStream.insert(byteStream.end(), data.begin(), data.end());

    return byteStream;
}

vector<int> unpack(vector<unsigned char> packedBytes, int numNeighbors)
{
    vector<int> unpackedValues;
    unpackedValues.reserve(numNeighbors);

    unpackedValues.push_back(int((packedBytes[0]) << 24 | (packedBytes[1]) << 16 | (packedBytes[2]) << 8 | (packedBytes[3])));

    int dataBytesLookedUp = 0;
    int numControlBytes = ceil(ceil((2.f * (float(numNeighbors) - 1.f)) / 8.f));
    for (size_t i = 0; i < numControlBytes; i++)
    {
        // Extract each 2-bit value from the byte
        for (int shiftAmount = 6; shiftAmount >= 0; shiftAmount -= 2)
        {
            unsigned char value = (packedBytes[i + 4] >> shiftAmount) & 0b11;  // Mask the last 2 bits

            // Convert the 2-bit value back to the original input range
            if (value == 0b00) {
                unpackedValues.push_back(0);
            }
            else if (value == 0b01) {
                unpackedValues.push_back(1);
            }
            else if (value == 0b10) {
                unpackedValues.push_back(int(packedBytes[dataBytesLookedUp + numControlBytes + 4]));
                dataBytesLookedUp++;
            }
            else if (value == 0b11) {
                int x = dataBytesLookedUp + numControlBytes + 4;
                unpackedValues.push_back(int((packedBytes[x]) << 24 | (packedBytes[x + 1]) << 16 | (packedBytes[x + 2]) << 8 | (packedBytes[x + 3])));
                dataBytesLookedUp += 4;
            }

            if (unpackedValues.size() >= numNeighbors) {
                break;  // Stop once we've unpacked the expected number of values
            }
        }
    }

    for (size_t i = 1; i < numNeighbors; i++)
    {
        unpackedValues[i] += unpackedValues[i - 1] + 1;
    }

    return unpackedValues;
}

int getMax(Particle* arr, int size)
{
    int m = arr[0].cellIndex;
    for (int i = 1; i < size; i++)
    {
        if (arr[i].cellIndex > m)
        {
            m = arr[i].cellIndex;
        }
    }
    return m;
}

void countingSort(Particle* arr, int size, int exp)
{
    vector<Particle> output(size);
    int count[10] = { 0 };

    for (int i = 0; i < size; i++) {
        int digit = static_cast<int>(arr[i].cellIndex / exp) % 10;
        count[digit]++;
    }

    for (int i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }

    for (int i = size - 1; i >= 0; i--) {
        int digit = static_cast<int>(arr[i].cellIndex / exp) % 10;
        output[count[digit] - 1] = arr[i];
        count[digit]--;
    }

    for (int i = 0; i < size; i++) {
        arr[i] = output[i];
    }
}

void radixSort(Particle* arr, int size)
{
    int maxCellIndex = getMax(arr, size);

    for (int exp = 1; maxCellIndex / exp > 0; exp *= 10)
    {
        countingSort(arr, size, exp);
    }
}