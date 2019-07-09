import {
    detectUtil,
    groupRectangles,
    convertRgbaToGrayscale,
    rescaleImage,
    computeSat,
    computeSquaredSat,
    computeRsat,
    compileClassifier,
    computeCanny
} from "./utils"

export interface Classifier extends Float32Array {
    tilted?: boolean
}

export const detector = (width: number, height: number, scaleFactor: number, classifier: Classifier) => {
    const canvas = document.createElement("canvas")
    canvas.width = width
    canvas.height = height
    const context = canvas.getContext("2d")!
    const tilted = classifier.tilted
    const numScales = ~~(Math.log(Math.min(width / classifier[0], height / classifier[1])) / Math.log(scaleFactor))
    let scaledGray = new Uint32Array(width * height)
    const compiledClassifiers: Classifier[] = []
    let scale = 1
    for (let i = 0; i < numScales; ++i) {
        const scaledWidth = ~~(width / scale)
        compiledClassifiers[i] = compileClassifier(classifier, scaledWidth)
        scale *= scaleFactor
    }

    let gray: Uint32Array | undefined
    let sat: Uint32Array | undefined
    let ssat: Uint32Array | undefined
    let rsat: Uint32Array | undefined
    let canny: Uint32Array | undefined
    let cannySat: Uint32Array | undefined

    const detect = (
        image: HTMLImageElement | HTMLCanvasElement | HTMLVideoElement,
        group: number = 1,
        stepSize: number = 1,
        roi?: number[],
        c?: any
    ) => {
        const w = canvas.width
        const h = canvas.height

        if (roi) context.drawImage(image, roi[0], roi[1], roi[2], roi[3], 0, 0, w, h)
        else context.drawImage(image, 0, 0, w, h)

        const imageData = context.getImageData(0, 0, w, h).data
        gray = convertRgbaToGrayscale(imageData)

        let rects: number[][] = []

        let s = 1
        for (let i = 0; i < numScales; ++i) {
            const scaledWidth = ~~(w / s)
            const scaledHeight = ~~(h / s)

            if (s === 1) {
                scaledGray.set(gray)
            } else {
                scaledGray = rescaleImage(gray, w, h, s, scaledGray)
            }

            if (c) {
                canny = computeCanny(scaledGray, scaledWidth, scaledHeight, canny)
                cannySat = computeSat(canny, scaledWidth, scaledHeight, cannySat)
            }

            sat = computeSat(scaledGray, scaledWidth, scaledHeight, sat)
            ssat = computeSquaredSat(scaledGray, scaledWidth, scaledHeight, ssat)
            if (tilted) rsat = computeRsat(scaledGray, scaledWidth, scaledHeight, rsat)

            const newRects = detectUtil(
                sat,
                rsat!,
                ssat,
                cannySat!,
                scaledWidth,
                scaledHeight,
                stepSize,
                compiledClassifiers[i]
            )
            for (let j = newRects.length - 1; j >= 0; --j) {
                newRects[j][0] *= s
                newRects[j][1] *= s
                newRects[j][2] *= s
                newRects[j][3] *= s
            }
            rects = rects.concat(newRects)

            s *= scaleFactor
        }
        return (group ? groupRectangles(rects, group) : rects).sort((r1, r2) => r2[4] - r1[4])
    }

    return { detect }
}
import { eye } from "./classifiers/eye"
import { frontalface } from "./classifiers/frontalface"

export const classifiers = { eye, frontalface }

export * from "./utils"
