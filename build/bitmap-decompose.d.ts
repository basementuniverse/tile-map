import { vec2 } from '@basementuniverse/vec';
export type Rectangle = {
    position: vec2;
    size: vec2;
};
export declare function bitmapToRectangles(bitmap: boolean[][]): Rectangle[];
