import { vec } from '@basementuniverse/vec';
export type Rectangle = {
    position: vec;
    size: vec;
};
export declare function bitmapToRectangles(bitmap: boolean[][]): Rectangle[];
