import { vec2 } from '@basementuniverse/vec';
import { Rectangle } from './bitmap-decompose';
export type TileMapOptionsData<T extends object = any> = Partial<Omit<TileMapOptions, 'preRender' | 'postRender' | 'preGenerateChunk' | 'postGenerateChunk' | 'debug'>> & {
    layers?: TileMapLayerOptionsData<T>[];
};
export type TileMapLayerOptionsData<T extends object = any> = Omit<TileMapLayerOptions, 'preRenderTile' | 'postRenderTile'> & {
    tiles?: (Omit<TileDefinition<T>, 'image'> & {
        imageName: string;
    })[];
    width?: number;
    data: number[];
};
export type TileMapOptions<T extends object = any> = {
    /**
     * The bounds of the tile map, measured in tiles
     *
     * Defines the position of the top-left and bottom-right corners of the
     * tile grid, used when rendering layer data
     *
     * If not defined, layer data will start at (0, 0)
     */
    bounds?: Bounds;
    /**
     * If true, the camera position will be clamped to the bounds
     *
     * Set this to false for infinite tilemaps
     *
     * Ignored if no bounds are defined
     *
     * Default is false
     */
    clampPositionToBounds: boolean;
    /**
     * The minimum scale factor
     */
    minScale?: number;
    /**
     * The maximum scale factor
     */
    maxScale?: number;
    /**
     * The size of each tile, measured in pixels
     *
     * Default is 16
     */
    tileSize: number;
    /**
     * A list of layers
     *
     * Defined in ascending render order; layers[0] will be rendered first, then
     * layers[1] on top of that, etc.
     */
    layers: TileMapLayerOptions<T>[];
    /**
     * The size of each render chunk, measured in tiles
     *
     * Default is 8
     */
    chunkSize: number;
    /**
     * Buffer area around the render area where we will load and render chunks,
     * measured in chunks
     *
     * This can be useful if rendering a chunk is quite slow; we can improve
     * the chances that a chunk will be ready to render by the time it appears
     * on-screen by increasing this number (although it means more chunks will
     * be rendered per frame)
     *
     * If set to a negative number, only render chunks which are fully inside
     * the screen bounds
     *
     * Default is 1
     */
    chunkBorder: number;
    /**
     * The maximum size of the LRU cache/queue
     *
     * Default is 64
     */
    chunkBufferMaxSize: number;
    /**
     * Optional hook called before rendering the tilemap
     *
     * @param context The context that the tilemap is being rendered into
     */
    preRender?: (context: CanvasRenderingContext2D, tileMap: TileMap<T>, screen: vec2, position: vec2, scale?: number) => void;
    /**
     * Optional hook called after rendering the tilemap
     *
     * @param context The context that the tilemap is being rendered into
     */
    postRender?: (context: CanvasRenderingContext2D, tileMap: TileMap<T>, screen: vec2, position: vec2, scale?: number) => void;
    /**
     * Optional hook called before generating a chunk
     *
     * Returns either a chunk, or a tuple containing a chunk and a boolean
     *
     * If the boolean element of the tuple is false, we skip default generation
     * even if tilemap data and images/tile definitions are defined
     *
     * In this case, the post-generate step will also be skipped
     *
     * @param context A context for the chunk's canvas
     */
    preGenerateChunk?: (context: CanvasRenderingContext2D, tileMap: TileMap<T>, tileBounds: Bounds, chunkPosition: vec2) => TileMapChunk | [TileMapChunk, boolean];
    /**
     * Optional hook called after generating a chunk
     *
     * @param canvas The chunk's canvas
     * @param context A context for the chunk's canvas
     */
    postGenerateChunk?: (canvas: HTMLCanvasElement, context: CanvasRenderingContext2D, tileMap: TileMap<T>, tileBounds: Bounds, chunkPosition: vec2) => TileMapChunk;
    /**
     * Optional debug options
     *
     * Can be a boolean value (in which case all sub-options will be set to the
     * same value), or an object allowing specific debug options to be enabled
     * individually
     */
    debug?: Partial<TileMapDebugOptions> | boolean;
};
export type TileMapLayerOptions<T extends object = any> = {
    /**
     * The name of this layer
     */
    name: string;
    /**
     * A list of tile definitions to use for tiles
     *
     * If this is not defined or empty, no tiles will be rendered
     *
     * The layer data will reference indexes in this array
     */
    tiles?: TileDefinition<T>[];
    /**
     * Layer data, represented as a 2d-array of indexes into the images array
     *
     * -1 means there is no tile at this position
     */
    data?: number[][];
    /**
     * Opacity of this layer, represented as a number between 0 (fully
     * transparent) and 1 (fully opaque)
     *
     * Default is 1
     */
    opacity?: number;
    /**
     * If true, each tile will be clipped to the tile size
     *
     * Default is false
     */
    clip?: boolean;
    /**
     * How should each tile's image be aligned within the tile?
     *
     * Default is TileAlignment.Center
     */
    alignment?: TileAlignment;
    /**
     * Optional hook called before rendering a tile
     *
     * @param context A context for the current chunk's canvas
     */
    preRenderTile?: (context: CanvasRenderingContext2D, tileMap: TileMap<T>, layer: TileMapLayerOptions<T>, chunkPosition: vec2, tilePosition: vec2) => void;
    /**
     * Optional hook called after rendering a tile
     *
     * If the tile doesn't have any data or a tile definition, the post-render
     * step will be skipped
     *
     * @param canvas The current chunk's canvas
     * @param context A context for the current chunk's canvas
     */
    postRenderTile?: (canvas: HTMLCanvasElement, context: CanvasRenderingContext2D, tileMap: TileMap<T>, layer: TileMapLayerOptions<T>, chunkPosition: vec2, tilePosition: vec2) => void;
};
type TileMapDebugOptions = {
    showOrigin: boolean;
    showChunkBorders: boolean;
    showChunkLabels: boolean;
    showTileBorders: boolean;
};
export type Bounds = {
    /**
     * The top-left corner of the tile map, measured in tiles
     */
    topLeft: vec2;
    /**
     * The bottom-right corner of the tile map, measured in tiles
     */
    bottomRight: vec2;
};
export declare enum TileAlignment {
    TopLeft = "top-left",
    Top = "top",
    TopRight = "top-right",
    Left = "left",
    Center = "center",
    Right = "right",
    BottomLeft = "bottom-left",
    Bottom = "bottom",
    BottomRight = "bottom-right"
}
export type TileDefinition<T extends object = any> = {
    name: string;
    image: TileMapImage;
    [key: string]: any;
} & T;
type TileMapImage = HTMLImageElement | HTMLCanvasElement;
type TileMapChunk = {
    chunkPosition: vec2;
    image: HTMLCanvasElement;
};
/**
 * Simplified interface for the camera component
 *
 * We can optionally pass this to the draw method instead of explicitly
 * passing the screen size, camera position and camera scale
 *
 * @see https://www.npmjs.com/package/@basementuniverse/camera
 */
interface Camera {
    position: vec2;
    readonly actualPosition: vec2;
    scale: number;
    readonly actualScale: number;
    bounds: {
        top: number;
        bottom: number;
        left: number;
        right: number;
    };
}
export declare function cameraBoundsToTileMapBounds(bounds: Camera['bounds']): Bounds;
export declare function cameraBoundsSize(bounds: Camera['bounds']): vec2;
export declare class TileMap<T extends object = any> {
    private static readonly DEFAULT_OPTIONS;
    private static readonly DEFAULT_LAYER_OPTIONS;
    private static readonly DEBUG_ORIGIN_COLOUR;
    private static readonly DEBUG_ORIGIN_LINE_WIDTH;
    private static readonly DEBUG_ORIGIN_SIZE;
    private static readonly DEBUG_CHUNK_BORDER_COLOUR;
    private static readonly DEBUG_CHUNK_BORDER_LINE_WIDTH;
    private static readonly DEBUG_CHUNK_LABEL_COLOUR;
    private static readonly DEBUG_CHUNK_LABEL_FONT;
    private static readonly DEBUG_TILE_BORDER_COLOUR;
    private static readonly DEBUG_TILE_BORDER_LINE_WIDTH;
    private options;
    private chunkBuffer;
    constructor(options?: Partial<TileMapOptions<T>>);
    /**
     * Get a (roughly minimal) set of rectangles which cover the tiles in a
     * given layer
     *
     * @param layerName The name of the layer to get rectangles for
     * @param fieldName We will check the truthyness of this field in the
     * tile definition
     * @param tileBounds Optional bounds to check within, relative to bounds
     * defined in options if any exist, otherwise relative to (0, 0)
     */
    getLayerRectangles(layerName: string, fieldName?: keyof TileDefinition<T>, tileBounds?: Bounds): Rectangle[];
    /**
     * Get the tile at a given position and in the specified layer
     *
     * If no layer is specified, return a dictionary of layer names to tile
     * definitions (i.e. return all layers)
     *
     * If no tile exists at this position, return null
     */
    getTileAtPosition(position: vec2, layerName?: string): TileDefinition<T> | null | {
        [name: string]: TileDefinition<T> | null;
    };
    private getTileAtPositionInLayer;
    private hashVector;
    draw(context: CanvasRenderingContext2D, camera: Camera): void;
    draw(context: CanvasRenderingContext2D, screen: vec2, position: vec2, scale: number): void;
    private performDraw;
    private generateChunk;
    private drawLine;
    private drawCross;
}
/**
 * Content Manager Processor wrapper which converts TileMapOptionsData into
 * TileMapOptions
 *
 * @see https://www.npmjs.com/package/@basementuniverse/content-manager
 */
export declare function tileMapOptionsContentProcessor<T extends object = any>(content: Record<string, {
    name: string;
    type: string;
    content: any;
    status: string;
}>, data: {
    name: string;
    type: string;
    content: any;
    status: string;
}, options?: Partial<{
    decompressData: boolean;
}>): Promise<void>;
export {};
