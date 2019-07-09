import {
    detectUtil,
    groupRectangles,
    convertRgbaToGrayscale,
    rescaleImage,
    computeSat,
    computeSquaredSat,
    computeRsat,
    compileClassifier
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

    const detect = (
        image: HTMLImageElement | HTMLCanvasElement | HTMLVideoElement,
        group: number,
        stepSize: number = 1
    ) => {
        const w = canvas.width
        const h = canvas.height

        context.drawImage(image, 0, 0, w, h)
        const imageData = context.getImageData(0, 0, w, h).data
        const gray = convertRgbaToGrayscale(imageData)
        let sat: Uint32Array | undefined
        let ssat: Uint32Array | undefined
        let rsat: Uint32Array | undefined
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

            sat = computeSat(scaledGray, scaledWidth, scaledHeight, sat)
            ssat = computeSquaredSat(scaledGray, scaledWidth, scaledHeight, ssat)
            if (tilted) rsat = computeRsat(scaledGray, scaledWidth, scaledHeight, rsat)

            const newRects = detectUtil(sat, rsat!, ssat, scaledWidth, scaledHeight, stepSize, compiledClassifiers[i])
            for (let j = newRects.length - 1; j >= 0; --j) {
                newRects[j][0] *= w / scaledWidth
                newRects[j][1] *= h / scaledHeight
                newRects[j][2] *= w / scaledWidth
                newRects[j][3] *= h / scaledHeight
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
