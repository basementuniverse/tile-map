# Game Component: Tile Map

Render 2d tile maps with layers, chunk rendering, tilesets, procedural generation, etc.

When this component is rendered (and at the same time given a position, screen size and scale), it will figure out which chunks should currently be visible on the screen. Then it will fetch these chunks from a cache, or generate them if they don't already exist. Then it renders them onto the screen.

By default we can provide tile definitions (each one has an image and optional extra data) and each layer can have tile data (a 2d array of indexes into the tile definitions array). This means that we can very quickly set up a simple layered tile map.

However, we also have various hooks (`preRender`, `preGenerateChunk`, `preRenderTile`, etc.) which allow us to customise the rendering process.

This makes it easy to do things like procedural generation, or custom post-processing effects like tile blending or marching squares, etc.

## Installation

```bash
npm install @basementuniverse/tile-map
```

## How to use

Create a tilemap:

```ts
import { TileMap } from '@basementuniverse/tile-map';

const tileMap = new TileMap({
  // options here...
})
```

Fetch the tile (or a map of tiles by layer name) at a position:

```ts
const tileAtPositionInLayer = tileMap.getTileAtPosition(
  { x: 0, y: 0 },
  'layer-name'
);
// returns a TileDefinition object or null

const tilesAtPosition = tileMap.getTilesAtPosition(
  { x: 0, y: 0 }
);
// returns an object where the keys are layer names and the values are TileDefinition objects or null
```

Render the tilemap every frame:

```ts
tileMap.draw(
  context,

  // a vec object representing the size of the viewport
  screen,

  // a vec object representing the camera position
  position,

  // a number representing the camera zoom level
  scale
);
```

_(see [here](https://www.npmjs.com/package/@basementuniverse/vec) for `vec` library)_

We can also fetch a list of rectangles representing the presence of a particular `TileDefinition` property in a given area (or across the whole tilemap).

This could be useful for collision detection, where checking collisions against potentially thousands of individual tiles could be expensive.

The set of rectangles returned isn't guaranteed to be the minimal set, but it's better than nothing :wink:

```ts
const collisionRectangles = tileMap.getLayerRectangles(
  layerName,

  // check this field in the TileDefinition, if it's truthy then a tile will form part of a rectangle
  fieldName,

  // optional, choose an area to fetch rectangles for
  tileBounds
);
// returns an array of rectangles

type Bounds = {
  topLeft: vec;
  bottomRight: vec;
};

type Rectangle = {
  position: vec;
  size: vec;
};
```

If you're using [this](https://www.npmjs.com/package/@basementuniverse/camera) camera component (or one with the same public interface), you can pass this to the `draw` method instead:

```ts
tileMap.draw(context, camera);
```

_Note: this will skip doing translate/scale transforms, assuming that the camera has already done them. It will also get the viewport/screen size from the camera. So, make sure to call `camera.draw()` before calling `tileMap.draw()`._

## Options

```ts
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
  preRender?: (
    context: CanvasRenderingContext2D,
    tileMap: TileMap<T>,
    screen: vec,
    position: vec,
    scale?: number
  ) => void;

  /**
   * Optional hook called after rendering the tilemap
   *
   * @param context The context that the tilemap is being rendered into
   */
  postRender?: (
    context: CanvasRenderingContext2D,
    tileMap: TileMap<T>,
    screen: vec,
    position: vec,
    scale?: number
  ) => void;

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
  preGenerateChunk?: (
    context: CanvasRenderingContext2D,
    tileMap: TileMap<T>,
    tileBounds: Bounds,
    chunkPosition: vec
  ) => TileMapChunk | [TileMapChunk, boolean];

  /**
   * Optional hook called after generating a chunk
   *
   * @param canvas The chunk's canvas
   * @param context A context for the chunk's canvas
   */
  postGenerateChunk?: (
    canvas: HTMLCanvasElement,
    context: CanvasRenderingContext2D,
    tileMap: TileMap<T>,
    tileBounds: Bounds,
    chunkPosition: vec
  ) => TileMapChunk;

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
  preRenderTile?: (
    context: CanvasRenderingContext2D,
    tileMap: TileMap<T>,
    layer: TileMapLayerOptions<T>,
    chunkPosition: vec,
    tilePosition: vec
  ) => void;

  /**
   * Optional hook called after rendering a tile
   *
   * If the tile doesn't have any data or a tile definition, the post-render
   * step will be skipped
   *
   * @param canvas The current chunk's canvas
   * @param context A context for the current chunk's canvas
   */
  postRenderTile?: (
    canvas: HTMLCanvasElement,
    context: CanvasRenderingContext2D,
    tileMap: TileMap<T>,
    layer: TileMapLayerOptions<T>,
    chunkPosition: vec,
    tilePosition: vec
  ) => void;
};
```

_(see `build/index.d.ts` for more details)_

## Content processor

A content processor function is provided for use with the [Content Manager](https://www.npmjs.com/package/@basementuniverse/content-manager).

```ts
import { tileMapOptionsContentProcessor } from '@basementuniverse/tile-map';
```

This function will take a JSON object resembling the `TileMapOptions` type (see below for an explanation of differences) and return a `TileMapOptions` object which can be passed into the `TileMap` constructor.

* Any `imageName: string` fields inside tile definitions in each layer will be replaced with an `image` field containing the image. The images will be fetched from the content manager.

* If we set the `decompressData` option, make sure there's a `width` field in each layer, and also make each layer's data field a 1d array of numbers, then the data will be decompressed using RLE and then split into rows.

```ts
ContentManager.initialise({
  processors: {
    tileMap: tileMapOptionsContentProcessor,
  },
});

ContentManager.load([
  {
    name: 'tile-map-1',
    type: 'json',
    args: [{
      // ...
      layers: [
        {
          // ...
          width: 4,
          data: [5, 0, 2, 1, 2, 0, 2, 1, 5, 0],
        }
      ]
    }],
    processors: [
      {
        name: 'tileMap',
        args: [{
          decompressData: true,
        }],
      },
    ],
  },
]);

// The layer data will be decompressed into:
// [
//   [0, 0, 0, 0],
//   [0, 1, 1, 0],
//   [0, 1, 1, 0],
//   [0, 0, 0, 0]
// ]
```

## Encoding tilemap data

A utility script is provided for RLE-encoding tilemap data.

```bash
node encode-rle "[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,3,3,0,0,0,0,0,0],[0,1,0,0,2,2,2,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,2,2,0,0,0,2,2,2,0,0,2],[0,0,1,1,0,2,0,0,0,2,0,0,0,2,2,0],[0,3,1,1,0,0,0,0,2,0,0,0,0,0,0,0],[0,0,0,1,1,0,0,0,2,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0],[0,0,1,1,0,0,0,2,0,0,0,0,0,0,0,0],[0,0,0,1,1,0,1,2,0,0,0,3,0,0,0,0],[0,0,0,0,0,1,0,0,2,0,0,0,0,0,0,0],[0,0,0,1,1,0,1,2,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,2,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,2,1,1,0,1,0,0,0,0,3,0],[0,0,0,0,2,0,0,0,1,3,3,0,0,0,0,0],[0,0,0,0,2,0,0,0,0,0,3,0,0,0,0,0]]"
```

This script will:

1. parse the data and flatten it into a 1d array
2. encode that array using RLE
3. dump the result to the console

```bash
[17,0,1,1,6,0,2,3,7,0,1,1,2,0,3,2,11,0,1,1,2,0,2,2,3,0,3,2,2,0,1,2,2,0,2,1,1,0,1,2,3,0,1,2,3,0,2,2,2,0,1,3,2,1,4,0,1,2,10,0,2,1,3,0,1,2,10,0,1,1,3,0,1,2,5,0,1,3,4,0,2,1,3,0,1,2,11,0,2,1,1,0,1,1,1,2,3,0,1,3,9,0,1,1,2,0,1,2,10,0,2,1,1,0,1,1,1,2,13,0,1,1,1,2,1,0,1,1,12,0,1,2,2,1,1,0,1,1,4,0,1,3,5,0,1,2,3,0,1,1,2,3,9,0,1,2,5,0,1,3,5,0]
```
