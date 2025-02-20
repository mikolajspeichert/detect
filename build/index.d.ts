export interface Classifier extends Float32Array {
    tilted?: boolean;
}
export declare const detector: (width: number, height: number, scaleFactor: number, classifier: Classifier) => {
    detect: (image: HTMLCanvasElement | HTMLImageElement | HTMLVideoElement, group?: number, stepSize?: number, roi?: number[] | undefined, c?: any) => number[][];
    canvas: HTMLCanvasElement;
};
export declare const classifiers: {
    eye: Classifier;
    frontalface: Classifier;
};
