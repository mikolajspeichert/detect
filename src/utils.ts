import { Classifier } from "."

export const convertRgbaToGrayscale = (src: Uint8ClampedArray, dst?: Uint32Array) => {
    const res = dst || new Uint32Array(src.length >> 2)

    for (let i = 0; i < src.length; i += 2) {
        res[i >> 2] = (src[i] * 4899 + src[++i] * 9617 + src[++i] * 1868 + 8192) >> 14
    }
    return res
}

export const rescaleImage = (
    src: Uint32Array,
    srcWidth: number,
    srcHeight: number,
    factor: number,
    dst?: Uint32Array
) => {
    const srcLength = srcHeight * srcWidth
    const dstWidth = ~~(srcWidth / factor)
    const dstHeight = ~~(srcHeight / factor)

    if (!dst) dst = new Uint32Array(dstWidth * srcHeight)

    for (let x = 0; x < dstWidth; ++x) {
        let di = x
        for (let srcIndex = ~~(x * factor), srcEnd = srcIndex + srcLength; srcIndex < srcEnd; srcIndex += srcWidth) {
            dst[di] = src[srcIndex]
            di += dstWidth
        }
    }

    let dstIndex = 0
    for (let y = 0, yEnd = dstHeight * factor; y < yEnd; y += factor) {
        for (let srcIndex = ~~y * dstWidth, srcEnd = srcIndex + dstWidth; srcIndex < srcEnd; ++srcIndex) {
            dst[dstIndex] = dst[srcIndex]
            ++dstIndex
        }
    }
    return dst
}

export const mirrorImage = (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array) => {
    if (!dst) dst = new Uint32Array(srcWidth * srcHeight)

    let index = 0
    for (let y = 0; y < srcHeight; ++y) {
        for (let x = srcWidth >> 1; x >= 0; --x) {
            const swap = src[index + x]
            dst[index + x] = src[index + srcWidth - 1 - x]
            dst[index + srcWidth - 1 - x] = swap
        }
        index += srcWidth
    }
    return dst
}

export const computeCanny = (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array) => {
    const srcLength = srcWidth * srcHeight
    if (!dst) dst = new Uint32Array(srcLength)
    const buffer1 = dst === src ? new Uint32Array(srcLength) : dst
    const buffer2 = new Uint32Array(srcLength)

    // Gaussian filter with size=5, sigma=sqrt(2) horizontal pass:
    for (let x = 2; x < srcWidth - 2; ++x) {
        let index = x
        for (let y = 0; y < srcHeight; ++y) {
            buffer1[index] =
                0.1117 * src[index - 2] +
                0.2365 * src[index - 1] +
                0.3036 * src[index] +
                0.2365 * src[index + 1] +
                0.1117 * src[index + 2]
            index += srcWidth
        }
    }

    // Gaussian filter with size=5, sigma=sqrt(2) vertical pass:
    for (let x = 0; x < srcWidth; ++x) {
        let index = x + srcWidth
        for (let y = 2; y < srcHeight - 2; ++y) {
            index += srcWidth
            buffer2[index] =
                0.1117 * buffer1[index - srcWidth - srcWidth] +
                0.2365 * buffer1[index - srcWidth] +
                0.3036 * buffer1[index] +
                0.2365 * buffer1[index + srcWidth] +
                0.1117 * buffer1[index + srcWidth + srcWidth]
        }
    }

    // Compute gradient:
    const abs = Math.abs
    for (let x = 2; x < srcWidth - 2; ++x) {
        let index = x + srcWidth
        for (let y = 2; y < srcHeight - 2; ++y) {
            index += srcWidth

            dst[index] =
                abs(
                    -buffer2[index - 1 - srcWidth] +
                        buffer2[index + 1 - srcWidth] -
                        2 * buffer2[index - 1] +
                        2 * buffer2[index + 1] -
                        buffer2[index - 1 + srcWidth] +
                        buffer2[index + 1 + srcWidth]
                ) +
                abs(
                    buffer2[index - 1 - srcWidth] +
                        2 * buffer2[index - srcWidth] +
                        buffer2[index + 1 - srcWidth] -
                        buffer2[index - 1 + srcWidth] -
                        2 * buffer2[index + srcWidth] -
                        buffer2[index + 1 + srcWidth]
                )
        }
    }
    return dst
}

export const computeSat = (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array) => {
    const dstWidth = srcWidth + 1

    if (!dst) dst = new Uint32Array(srcWidth * srcHeight + dstWidth + srcHeight)

    for (let i = srcHeight * dstWidth; i >= 0; i -= dstWidth) dst[i] = 0

    for (let x = 1; x <= srcWidth; ++x) {
        let column_sum = 0
        let index = x
        dst[x] = 0

        for (let y = 1; y <= srcHeight; ++y) {
            column_sum += src[index - y]
            index += dstWidth
            dst[index] = dst[index - 1] + column_sum
        }
    }
    return dst
}

export const computeSquaredSat = (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array) => {
    const dstWidth = srcWidth + 1

    if (!dst) dst = new Uint32Array(srcWidth * srcHeight + dstWidth + srcHeight)

    for (let i = srcHeight * dstWidth; i >= 0; i -= dstWidth) dst[i] = 0

    for (let x = 1; x <= srcWidth; ++x) {
        let column_sum = 0
        let index = x
        dst[x] = 0
        for (let y = 1; y <= srcHeight; ++y) {
            const val = src[index - y]
            column_sum += val * val
            index += dstWidth
            dst[index] = dst[index - 1] + column_sum
        }
    }
    return dst
}

export const computeRsat = (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array) => {
    const dstWidth = srcWidth + 1
    const srcHeightTimesDstWidth = srcHeight * dstWidth

    if (!dst) dst = new Uint32Array(srcWidth * srcHeight + dstWidth + srcHeight)

    for (let i = srcHeightTimesDstWidth; i >= 0; i -= dstWidth) dst[i] = 0

    for (let i = 0; i < dstWidth; ++i) dst[i] = 0

    let index = 0
    for (let y = 0; y < srcHeight; ++y) {
        for (let x = 0; x < srcWidth; ++x) {
            dst[index + dstWidth + 1] = src[index - y] + dst[index]
            ++index
        }
        dst[index + dstWidth] += dst[index]
        index++
    }

    for (let x = srcWidth - 1; x > 0; --x) {
        let i = x + srcHeightTimesDstWidth
        for (let y = srcHeight; y > 0; --y) {
            i -= dstWidth
            dst[i + dstWidth] += dst[i] + dst[i + 1]
        }
    }

    return dst
}

export const equalizeHistogram = (src: number[], step: number = 5, dst = src) => {
    const srcLength = src.length

    // Compute histogram and histogram sum:
    const hist = Array(256).fill(0)

    for (let i = 0; i < srcLength; i += step) {
        ++hist[src[i]]
    }

    // Compute integral histogram:
    const norm = (255 * step) / srcLength
    let prev = 0
    for (let i = 0; i < 256; ++i) {
        let h = hist[i]
        prev = h += prev
        hist[i] = h * norm // For non-integer src: ~~(h * norm + 0.5);
    }

    // Equalize image:
    for (let i = 0; i < srcLength; ++i) {
        dst[i] = hist[src[i]]
    }
    return dst
}

export const mirrorClassifier = (src: Uint32Array, dst: Uint32Array) => {
    if (!dst) dst = new Uint32Array(src)
    const windowWidth = src[0]

    for (let i = 1, iEnd = src.length - 1; i < iEnd; ) {
        ++i
        for (let j = 0, jEnd = src[++i]; j < jEnd; ++j) {
            if (src[++i]) {
                // Simple classifier is tilted:
                for (const kEnd = i + src[++i] * 5; i < kEnd; ) {
                    dst[i + 1] = windowWidth - src[i + 1]
                    const width = src[i + 3]
                    dst[i + 3] = src[i + 4]
                    dst[i + 4] = width
                    i += 5
                }
            } else {
                // Simple classifier is not tilted:
                for (const kEnd = i + src[++i] * 5; i < kEnd; ) {
                    dst[i + 1] = windowWidth - src[i + 1] - src[i + 3]
                    i += 5
                }
            }
            i += 3
        }
    }
    return dst
}

export const compileClassifier = (src: Classifier, width: number, dst?: Classifier) => {
    width += 1
    if (!dst) dst = new Float32Array(src.length)
    const dstUint32 = new Uint32Array(dst.buffer)

    dstUint32[0] = src[0]
    dstUint32[1] = src[1]
    let dstIndex = 1
    for (let srcIndex = 1, iEnd = src.length - 1; srcIndex < iEnd; ) {
        dst[++dstIndex] = src[++srcIndex]

        const numComplexClassifiers = (dstUint32[++dstIndex] = src[++srcIndex])
        for (let j = 0, jEnd = numComplexClassifiers; j < jEnd; ++j) {
            const tilted = (dst[++dstIndex] = src[++srcIndex])
            const numFeaturesTimes3 = (dstUint32[++dstIndex] = src[++srcIndex] * 3)
            if (tilted) {
                for (const kEnd = dstIndex + numFeaturesTimes3; dstIndex < kEnd; ) {
                    dstUint32[++dstIndex] = src[++srcIndex] + src[++srcIndex] * width
                    dstUint32[++dstIndex] = src[++srcIndex] * (width + 1) + ((src[++srcIndex] * (width - 1)) << 16)
                    dst[++dstIndex] = src[++srcIndex]
                }
            } else {
                for (const kEnd = dstIndex + numFeaturesTimes3; dstIndex < kEnd; ) {
                    dstUint32[++dstIndex] = src[++srcIndex] + src[++srcIndex] * width
                    dstUint32[++dstIndex] = src[++srcIndex] + ((src[++srcIndex] * width) << 16)
                    dst[++dstIndex] = src[++srcIndex]
                }
            }

            const inverseClassifierThreshold = 1 / src[++srcIndex]
            for (let k = 0; k < numFeaturesTimes3; ) {
                dst[dstIndex - k] *= inverseClassifierThreshold
                k += 3
            }

            if (inverseClassifierThreshold < 0) {
                dst[dstIndex + 2] = src[++srcIndex]
                dst[dstIndex + 1] = src[++srcIndex]
                dstIndex += 2
            } else {
                dst[++dstIndex] = src[++srcIndex]
                dst[++dstIndex] = src[++srcIndex]
            }
        }
    }
    return dst.subarray(0, dstIndex + 1)
}

export const detectUtil = (
    sat: Uint32Array,
    rsat: Uint32Array,
    ssat: Uint32Array,
    cannySat: Uint32Array,
    width: number,
    height: number,
    step: number,
    classifier: Float32Array
) => {
    width += 1
    height += 1

    const classifierUint32 = new Uint32Array(classifier.buffer)
    const windowWidth = classifierUint32[0]
    const windowHeight = classifierUint32[1]
    const windowHeightTimesWidth = windowHeight * width
    const area = windowWidth * windowHeight
    const inverseArea = 1 / area
    const widthTimesStep = width * step
    const rects = []

    for (let x = 0; x + windowWidth < width; x += step) {
        let satIndex = x
        for (let y = 0; y + windowHeight < height; y += step) {
            const satIndex1 = satIndex + windowWidth
            const satIndex2 = satIndex + windowHeightTimesWidth
            const satIndex3 = satIndex2 + windowWidth

            if (cannySat) {
                const edgesDensity =
                    (cannySat[satIndex] - cannySat[satIndex1] - cannySat[satIndex2] + cannySat[satIndex3]) * inverseArea
                if (edgesDensity < 60 || edgesDensity > 200) {
                    satIndex += widthTimesStep
                    continue
                }
            }

            // Normalize mean and variance of window area:
            const mean = sat[satIndex] - sat[satIndex1] - sat[satIndex2] + sat[satIndex3]
            const variance = (ssat[satIndex] - ssat[satIndex1] - ssat[satIndex2] + ssat[satIndex3]) * area - mean * mean
            const std = variance > 1 ? Math.sqrt(variance) : 1
            let found = true

            // Evaluate cascade classifier aka 'stages':
            for (let i = 1, iEnd = classifier.length - 1; i < iEnd; ) {
                const complexClassifierThreshold = classifier[++i]
                // Evaluate complex classifiers aka 'trees':
                let complexClassifierSum = 0
                for (let j = 0, jEnd = classifierUint32[++i]; j < jEnd; ++j) {
                    // Evaluate simple classifiers aka 'nodes':
                    let simpleClassifierSum = 0

                    if (classifierUint32[++i]) {
                        // Simple classifier is tilted:
                        for (const kEnd = i + classifierUint32[++i]; i < kEnd; ) {
                            const f1 = satIndex + classifierUint32[++i]
                            const packed = classifierUint32[++i]
                            const f2 = f1 + (packed & 0xffff)
                            const f3 = f1 + ((packed >> 16) & 0xffff)

                            simpleClassifierSum +=
                                classifier[++i] * (rsat[f1] - rsat[f2] - rsat[f3] + rsat[f2 + f3 - f1])
                        }
                    } else {
                        // Simple classifier is not tilted:
                        for (const kEnd = i + classifierUint32[++i]; i < kEnd; ) {
                            const f1 = satIndex + classifierUint32[++i]
                            const packed = classifierUint32[++i]
                            const f2 = f1 + (packed & 0xffff)
                            const f3 = f1 + ((packed >> 16) & 0xffff)

                            simpleClassifierSum += classifier[++i] * (sat[f1] - sat[f2] - sat[f3] + sat[f2 + f3 - f1])
                        }
                    }
                    complexClassifierSum += classifier[i + (simpleClassifierSum > std ? 2 : 1)]
                    i += 2
                }
                if (complexClassifierSum < complexClassifierThreshold) {
                    found = false
                    break
                }
            }
            if (found) rects.push([x, y, windowWidth, windowHeight])
            satIndex += widthTimesStep
        }
    }
    return rects
}

export const groupRectangles = (rects: number[][], minNeighbors: number, confluence: number = 0.25) => {
    const rectsLength = rects.length

    // Partition rects into similarity classes:
    let numClasses = 0
    const labels = new Array(rectsLength)
    for (let i = 0; i < rectsLength; ++i) {
        labels[i] = 0
    }

    for (let i = 0; i < rectsLength; ++i) {
        let found = false
        for (let j = 0; j < i; ++j) {
            // Determine similarity:
            const rect1 = rects[i]
            const rect2 = rects[j]
            const delta = confluence * (Math.min(rect1[2], rect2[2]) + Math.min(rect1[3], rect2[3]))
            if (
                Math.abs(rect1[0] - rect2[0]) <= delta &&
                Math.abs(rect1[1] - rect2[1]) <= delta &&
                Math.abs(rect1[0] + rect1[2] - rect2[0] - rect2[2]) <= delta &&
                Math.abs(rect1[1] + rect1[3] - rect2[1] - rect2[3]) <= delta
            ) {
                labels[i] = labels[j]
                found = true
                break
            }
        }
        if (!found) {
            labels[i] = numClasses++
        }
    }

    // Compute average rectangle (group) for each cluster:
    const groups: number[][] = new Array(numClasses)

    for (let i = 0; i < numClasses; ++i) {
        groups[i] = [0, 0, 0, 0, 0]
    }

    for (let i = 0; i < rectsLength; ++i) {
        const rect = rects[i]
        const group = groups[labels[i]]
        group[0] += rect[0]
        group[1] += rect[1]
        group[2] += rect[2]
        group[3] += rect[3]
        ++group[4]
    }

    for (let i = 0; i < numClasses; ++i) {
        let numNeighbors = groups[i][4]
        if (numNeighbors >= minNeighbors) {
            const group = groups[i]
            numNeighbors = 1 / numNeighbors
            group[0] *= numNeighbors
            group[1] *= numNeighbors
            group[2] *= numNeighbors
            group[3] *= numNeighbors
        } else groups.splice(i, 1)
    }

    // Filter out small rectangles inside larger rectangles:
    // Filter out small rectangles inside larger rectangles:
    const filteredGroups = []
    for (let i = 0; i < numClasses; ++i) {
        const r1 = groups[i]
        let j = 0
        for (j = i + 1; j < numClasses; ++j) {
            const r2 = groups[j]

            const dx = r2[2] * confluence // * 0.2;
            const dy = r2[3] * confluence // * 0.2;

            // Not antisymmetric, must check both r1 > r2 and r2 > r1:
            if (
                (r1[0] >= r2[0] - dx &&
                    r1[1] >= r2[1] - dy &&
                    r1[0] + r1[2] <= r2[0] + r2[2] + dx &&
                    r1[1] + r1[3] <= r2[1] + r2[3] + dy) ||
                (r2[0] >= r1[0] - dx &&
                    r2[1] >= r1[1] - dy &&
                    r2[0] + r2[2] <= r1[0] + r1[2] + dx &&
                    r2[1] + r2[3] <= r1[1] + r1[3] + dy)
            ) {
                break
            }
        }

        if (j === numClasses) {
            filteredGroups.push(r1)
        }
    }
    return filteredGroups
}

export const Smoother = (alphas: number[], initialValues: number[], lookAhead: number = 1.0) => {
    let lastUpdate = new Date().getTime()
    const initialAlphas = [...alphas]
    const as = [...alphas]
    const a = [...initialValues]
    const b = [...initialValues]
    const numValues = initialValues.length

    const smooth = (values: number[]) => {
        const smoothedValues = []

        // time in seconds since last update:
        let time = new Date().getTime() - lastUpdate
        lastUpdate += time
        time /= 1000

        // update:
        for (let i = 0; i < numValues; ++i) {
            // Wright's modification of Holt's method for irregular data:
            as[i] = as[i] / (as[i] + Math.pow(1 - initialAlphas[i], time))

            const oldA = a[i]
            a[i] = as[i] * values[i] + (1 - as[i]) * (a[i] + b[i] * time)
            b[i] = (as[i] * (a[i] - oldA)) / time + (1 - as[i]) * b[i]

            smoothedValues[i] = a[i] + time * lookAhead * b[i]

            // Alternative approach:
            // a[i] = as[i] * values[i] + (1 - as[i]) * a[i];
            // b[i] = as[i] * a[i] + (1 - as[i]) * b[i];
            // smoothedValues[i] = 2*a[i] - 1*b[i];*/
        }

        return smoothedValues
    }
    return { smooth }
}
