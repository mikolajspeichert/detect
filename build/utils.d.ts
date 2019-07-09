import { Classifier } from ".";
export declare const convertRgbaToGrayscale: (src: Uint8ClampedArray, dst?: Uint32Array | undefined) => Uint32Array;
export declare const rescaleImage: (src: Uint32Array, srcWidth: number, srcHeight: number, factor: number, dst?: Uint32Array | undefined) => Uint32Array;
export declare const mirrorImage: (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array | undefined) => Uint32Array;
export declare const computeCanny: (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array | undefined) => Uint32Array;
export declare const computeSat: (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array | undefined) => Uint32Array;
export declare const computeSquaredSat: (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array | undefined) => Uint32Array;
export declare const computeRsat: (src: Uint32Array, srcWidth: number, srcHeight: number, dst?: Uint32Array | undefined) => Uint32Array;
export declare const equalizeHistogram: (src: number[], step?: number, dst?: number[]) => number[];
export declare const mirrorClassifier: (src: Uint32Array, dst: Uint32Array) => Uint32Array;
export declare const compileClassifier: (src: Classifier, width: number, dst?: Classifier | undefined) => Float32Array;
export declare const detectUtil: (sat: Uint32Array, rsat: Uint32Array, ssat: Uint32Array, cannySat: Uint32Array, width: number, height: number, step: number, classifier: Float32Array) => number[][];
export declare const groupRectangles: (rects: number[][], minNeighbors: number, confluence?: number) => number[][];
export declare const Smoother: (alphas: number[], initialValues: number[], lookAhead?: number) => {
    smooth: (values: number[]) => number[];
};
