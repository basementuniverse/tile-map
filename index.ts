import { vec } from '@basementuniverse/vec';
import { LRUMap } from 'lru_map';

export type TileMapOptionsData = Partial<Omit<
  TileMapOptions,
  | 'preRender'
  | 'postRender'
  | 'preGenerateChunk'
  | 'postGenerateChunk'
  | 'debug'
>> & {
  layers?: TileMapLayerOptionsData[];
};

export type TileMapLayerOptionsData = Omit<
  TileMapLayerOptions,
  | 'preRenderTile'
  | 'postRenderTile'
> & {
  tiles?: Omit<TileDefinition, 'image'> & {
    imageName: string;
    [key: string]: any;
  }[];
};

export type TileMapOptions<
  T extends TileDefinition | undefined = undefined
> = {
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

export type TileMapLayerOptions<
  T extends TileDefinition | undefined = undefined
> = {
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
  tiles?: T[];

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
  topLeft: vec;

  /**
   * The bottom-right corner of the tile map, measured in tiles
   */
  bottomRight: vec;
};

export enum TileAlignment {
  TopLeft = 0,
  Top,
  TopRight,
  Left,
  Center,
  Right,
  BottomLeft,
  Bottom,
  BottomRight,
}

export type TileDefinition = {
  name: string;
  image: TileMapImage;
};

type TileMapImage = HTMLImageElement | HTMLCanvasElement;

type TileMapChunk = {
  chunkPosition: vec;
  image: HTMLCanvasElement;
};

function clamp(a: number, min: number = 0, max: number = 1): number {
  return a < min ? min : (a > max ? max : a);
}

function pointInRectangle(
  point: vec,
  topLeft: vec,
  bottomRight: vec
): boolean {
  return (
    point.x >= topLeft.x &&
    point.y >= topLeft.y &&
    point.x < bottomRight.x &&
    point.y < bottomRight.y
  );
}

export class TileMap<T extends TileDefinition | undefined = undefined> {
  // TODO

  // https://www.npmjs.com/package/rectangle-decomposition
  // https://www.npmjs.com/package/fast-rle

  private static readonly DEFAULT_OPTIONS: TileMapOptions = {
    clampPositionToBounds: true,
    tileSize: 16,
    layers: [
      {
        name: 'default',
      },
    ],
    chunkSize: 8,
    chunkBorder: 1,
    chunkBufferMaxSize: 64,
  };

  private static readonly DEBUG_ORIGIN_COLOUR = 'cyan';
  private static readonly DEBUG_ORIGIN_LINE_WIDTH = 2;
  private static readonly DEBUG_ORIGIN_SIZE = 10;

  private static readonly DEBUG_CHUNK_BORDER_COLOUR = 'yellow';
  private static readonly DEBUG_CHUNK_BORDER_LINE_WIDTH = 2;

  private static readonly DEBUG_CHUNK_LABEL_COLOUR = 'white';
  private static readonly DEBUG_CHUNK_LABEL_FONT = '12px monospace';

  private static readonly DEBUG_TILE_BORDER_COLOUR = 'orange';
  private static readonly DEBUG_TILE_BORDER_LINE_WIDTH = 1;

  private options: TileMapOptions<T> & {
    debug: Required<TileMapDebugOptions>;
  };

  private chunkBuffer: LRUMap<string, TileMapChunk>;

  public constructor(options?: Partial<TileMapOptions<T>>) {
    const actualOptions = Object.assign(
      {},
      TileMap.DEFAULT_OPTIONS,
      options ?? {}
    );

    if (!actualOptions.debug || actualOptions.debug === true) {
      actualOptions.debug = {
        showOrigin: !!actualOptions.debug,
        showChunkBorders: !!actualOptions.debug,
        showChunkLabels: !!actualOptions.debug,
        showTileBorders: !!actualOptions.debug,
      };
    }

    this.options = actualOptions as typeof this.options;

    this.chunkBuffer = new LRUMap(this.options.chunkBufferMaxSize);
  }

  /**
   * Get a minimal set of rectangles which cover the tiles in a given layer
   *
   * @param layerName The name of the layer to get rectangles for
   * @param fieldName We will check the truthyness of this field in the
   * tile definition
   * @param tileBounds Optional bounds to check
   */
  public getLayerRectangles(
    layerName: string,
    fieldName?: string,
    tileBounds?: Bounds
  ): any {
    // TODO
  }

  /**
   * Get the tile at a given position and in the specified layer
   *
   * If no layer is specified, return a dictionary of layer names to tile
   * definitions (i.e. return all layers)
   *
   * If no tile exists at this position, return null
   */
  public getTileAtPosition(
    position: vec,
    layerName?: string
  ): T | null | { [name: string]: T | null } {
    if (layerName) {
      return this.getTileAtPositionInLayer(position, layerName);
    }

    const result: { [name: string]: T | null } = {};
    for (const layer of this.options.layers) {
      result[layer.name] = this.getTileAtPositionInLayer(position, layer.name);
    }

    return result;
  }

  private getTileAtPositionInLayer(
    position: vec,
    layerName: string
  ): T | null {
    const tilePosition = vec.map(
      vec.mul(position, 1 / this.options.tileSize),
      Math.floor
    );

    const layer = this.options.layers.find((l) => l.name === layerName);
    if (!layer) {
      return null;
    }

    const tileData = layer.data?.[tilePosition.y]?.[tilePosition.x];
    if (tileData === undefined || tileData === -1) {
      return null;
    }

    if (layer.tiles) {
      return layer.tiles[tileData] ?? null;
    }

    return null;
  }

  private hashVector(v: vec): string {
    return vec.str(v);
  }

  public draw(
    context: CanvasRenderingContext2D,
    screen: vec,
    position: vec,
    scale: number
  ) {
    const absoluteChunkSize = this.options.tileSize * this.options.chunkSize;
    const chunkBorder = vec(this.options.chunkBorder);

    // Maybe clamp scale
    let actualScale = scale;
    if (this.options.minScale && actualScale < this.options.minScale) {
      actualScale = this.options.minScale;
    }
    if (this.options.maxScale && actualScale > this.options.maxScale) {
      actualScale = this.options.maxScale;
    }

    // Maybe clamp position to bounds
    let actualPosition = vec(position);
    if (this.options.bounds && this.options.clampPositionToBounds) {
      const tileSizeScaled = this.options.tileSize / actualScale;
      const halfScreenScaled = vec.map(
        vec.mul(screen, 1 / (actualScale * 2)),
        Math.ceil
      );
      const minPosition = vec(
        this.options.bounds.topLeft.x * tileSizeScaled + halfScreenScaled.x,
        this.options.bounds.topLeft.y * tileSizeScaled + halfScreenScaled.y
      );
      const maxPosition = vec(
        this.options.bounds.bottomRight.x * tileSizeScaled - halfScreenScaled.x,
        this.options.bounds.bottomRight.y * tileSizeScaled - halfScreenScaled.y
      );

      actualPosition = vec(
        clamp(actualPosition.x, minPosition.x, maxPosition.x),
        clamp(actualPosition.y, minPosition.y, maxPosition.y)
      );
    }

    const screenSizeInChunks = vec.map(
      vec.mul(
        screen,
        1 / (absoluteChunkSize * actualScale)
      ),
      Math.ceil
    );
    const screenCenterChunk = vec.map(
      vec.mul(actualPosition, 1 / absoluteChunkSize),
      Math.floor
    );
    const topLeftChunk = vec.sub(
      vec.sub(
        screenCenterChunk,
        vec.map(
          vec.mul(screenSizeInChunks, 0.5),
          Math.ceil
        )
      ),
      chunkBorder
    );
    const bottomRightChunk = vec.add(
      vec.add(
        screenCenterChunk,
        vec.map(
          vec.mul(screenSizeInChunks, 0.5),
          Math.ceil
        )
      ),
      chunkBorder
    );

    context.save();
    context.scale(actualScale, actualScale);
    context.translate(
      -actualPosition.x + screen.x / (actualScale * 2),
      -actualPosition.y + screen.y / (actualScale * 2)
    );

    this.options.preRender?.(
      context,
      this,
      screen,
      actualPosition,
      actualScale
    );

    // Render chunks
    for (let y = topLeftChunk.y; y < bottomRightChunk.y; y++) {
      for (let x = topLeftChunk.x; x < bottomRightChunk.x; x++) {
        const chunkPosition = vec(x, y);
        const chunkAbsolutePosition = vec.mul(chunkPosition, absoluteChunkSize);

        // Check if we have this chunk in the cache
        const chunkHash = this.hashVector(chunkPosition);
        if (!this.chunkBuffer.has(chunkHash)) {
          this.chunkBuffer.set(chunkHash, this.generateChunk(
            chunkPosition,
            absoluteChunkSize
          ));
        }

        const chunk = this.chunkBuffer.get(chunkHash);
        if (chunk) {
          context.drawImage(
            chunk.image,
            chunkAbsolutePosition.x,
            chunkAbsolutePosition.y
          );
        }
      }
    }

    this.options.postRender?.(
      context,
      this,
      screen,
      actualPosition,
      actualScale
    );

    // Render debug helpers
    if (this.options.debug.showTileBorders) {
      const topLeftTile = vec.mul(
        vec.sub(
          screenCenterChunk,
          vec.add(
            vec.map(
              vec.mul(screenSizeInChunks, 0.5),
              Math.ceil
            ),
            vec(1)
          )
        ),
        this.options.chunkSize
      );
      const bottomRightTile = vec.mul(
        vec.add(
          screenCenterChunk,
          vec.add(
            vec.map(
              vec.mul(screenSizeInChunks, 0.5),
              Math.ceil
            ),
            vec(1)
          )
        ),
        this.options.chunkSize
      );

      for (let y = topLeftTile.y; y < bottomRightTile.y; y++) {
        this.drawLine(
          context,
          vec(
            actualPosition.x - screen.x / (actualScale * 2),
            y * this.options.tileSize
          ),
          vec(
            actualPosition.x + screen.x / (actualScale * 2),
            y * this.options.tileSize
          ),
          TileMap.DEBUG_TILE_BORDER_COLOUR,
          TileMap.DEBUG_TILE_BORDER_LINE_WIDTH
        );
      }
      for (let x = topLeftTile.x; x < bottomRightTile.x; x++) {
        this.drawLine(
          context,
          vec(
            x * this.options.tileSize,
            actualPosition.y - screen.y / (actualScale * 2)
          ),
          vec(
            x * this.options.tileSize,
            actualPosition.y + screen.y / (actualScale * 2)
          ),
          TileMap.DEBUG_TILE_BORDER_COLOUR,
          TileMap.DEBUG_TILE_BORDER_LINE_WIDTH
        );
      }
    }

    if (this.options.debug.showChunkBorders) {
      for (let y = topLeftChunk.y; y < bottomRightChunk.y; y++) {
        this.drawLine(
          context,
          vec(
            actualPosition.x - screen.x / (actualScale * 2),
            y * absoluteChunkSize
          ),
          vec(
            actualPosition.x + screen.x / (actualScale * 2),
            y * absoluteChunkSize
          ),
          TileMap.DEBUG_CHUNK_BORDER_COLOUR,
          TileMap.DEBUG_CHUNK_BORDER_LINE_WIDTH
        );
      }
      for (let x = topLeftChunk.x; x < bottomRightChunk.x; x++) {
        this.drawLine(
          context,
          vec(
            x * absoluteChunkSize,
            actualPosition.y - screen.y / (actualScale * 2)
          ),
          vec(
            x * absoluteChunkSize,
            actualPosition.y + screen.y / (actualScale * 2)
          ),
          TileMap.DEBUG_CHUNK_BORDER_COLOUR,
          TileMap.DEBUG_CHUNK_BORDER_LINE_WIDTH
        );
      }
    }

    if (this.options.debug.showChunkLabels) {
      context.save();
      context.fillStyle = TileMap.DEBUG_CHUNK_LABEL_COLOUR;
      context.font = TileMap.DEBUG_CHUNK_LABEL_FONT;
      context.textBaseline = 'middle';
      context.textAlign = 'center';

      for (let y = topLeftChunk.y; y < bottomRightChunk.y; y++) {
        for (let x = topLeftChunk.x; x < bottomRightChunk.x; x++) {
          context.fillText(
            `${x}, ${y}`,
            x * absoluteChunkSize + absoluteChunkSize / 2,
            y * absoluteChunkSize + absoluteChunkSize / 2
          );
        }
      }

      context.restore();
    }

    if (
      this.options.debug.showOrigin &&
      pointInRectangle(vec(0, 0), topLeftChunk, bottomRightChunk)
    ) {
      this.drawCross(
        context,
        vec(0, 0),
        TileMap.DEBUG_ORIGIN_COLOUR,
        TileMap.DEBUG_ORIGIN_LINE_WIDTH,
        TileMap.DEBUG_ORIGIN_SIZE
      );
    }

    context.restore();
  }

  private generateChunk(
    chunkPosition: vec,
    absoluteChunkSize: number
  ): TileMapChunk {
    const chunkCanvas = document.createElement('canvas');
    const chunkContext = chunkCanvas.getContext('2d')!;

    chunkCanvas.width = absoluteChunkSize;
    chunkCanvas.height = absoluteChunkSize;

    let chunk: TileMapChunk = {
      chunkPosition,
      image: chunkCanvas,
    };

    const topLeftTile = vec.mul(chunkPosition, this.options.chunkSize);
    const bottomRightTile = vec.add(
      topLeftTile,
      vec(this.options.chunkSize - 1)
    );
    const boundsTopLeft = this.options.bounds?.topLeft ?? vec(0);

    if (this.options.preGenerateChunk) {
      const result = this.options.preGenerateChunk(
        chunkContext,
        this,
        {
          topLeft: topLeftTile,
          bottomRight: bottomRightTile,
        },
        chunkPosition
      );

      if (Array.isArray(result)) {
        if (!result[1]) {
          return chunk;
        }
      }
    }

    // Default generation, render tiles from tilemap data
    for (const layer of this.options.layers) {
      chunkContext.save();
      chunkContext.globalAlpha = layer.opacity ?? 1;

      const alignment = layer.alignment ?? TileAlignment.Center;

      for (let y = topLeftTile.y; y <= bottomRightTile.y; y++) {
        for (let x = topLeftTile.x; x <= bottomRightTile.x; x++) {
          const tilePosition = vec(x, y);

          layer.preRenderTile?.(
            chunkContext,
            this,
            layer,
            chunkPosition,
            tilePosition
          );

          const tileDataPosition = vec.sub(
            tilePosition,
            boundsTopLeft
          );

          if (tileDataPosition.x < 0 || tileDataPosition.y < 0) {
            continue;
          }

          const tileData = layer.data
            ?.[tileDataPosition.y]
            ?.[tileDataPosition.x];
          if (tileData === undefined || tileData === -1) {
            continue;
          }

          const tileImage = layer.tiles?.[tileData]?.image;
          if (!tileImage) {
            continue;
          }

          const tileAbsolutePosition = vec.sub(
            vec.mul(
              tilePosition,
              this.options.tileSize
            ),
            vec.mul(chunkPosition, absoluteChunkSize)
          );

          // Tile clipping
          if (layer.clip) {
            chunkContext.save();
            chunkContext.beginPath();
            chunkContext.rect(
              tileAbsolutePosition.x,
              tileAbsolutePosition.y,
              this.options.tileSize,
              this.options.tileSize
            );
            chunkContext.clip();
          }

          // Tile alignment
          let tileImageAbsolutePosition: vec;
          switch (alignment) {
            case TileAlignment.TopLeft:
              tileImageAbsolutePosition = vec(tileAbsolutePosition);
              break;

            case TileAlignment.Top:
              tileImageAbsolutePosition = vec(
                (
                  tileAbsolutePosition.x + this.options.tileSize / 2
                ) - tileImage.width / 2,
                tileAbsolutePosition.y
              );
              break;

            case TileAlignment.TopRight:
              tileImageAbsolutePosition = vec(
                tileAbsolutePosition.x + this.options.tileSize - tileImage.width,
                tileAbsolutePosition.y
              );
              break;

            case TileAlignment.Left:
              tileImageAbsolutePosition = vec(
                tileAbsolutePosition.x,
                (
                  tileAbsolutePosition.y + this.options.tileSize / 2
                ) - tileImage.height / 2
              );
              break;

            case TileAlignment.Center:
              tileImageAbsolutePosition = vec(
                (
                  tileAbsolutePosition.x + this.options.tileSize / 2
                ) - tileImage.width / 2,
                (
                  tileAbsolutePosition.y + this.options.tileSize / 2
                ) - tileImage.height / 2
              );
              break;

            case TileAlignment.Right:
              tileImageAbsolutePosition = vec(
                tileAbsolutePosition.x + this.options.tileSize - tileImage.width,
                (
                  tileAbsolutePosition.y + this.options.tileSize / 2
                ) - tileImage.height / 2
              );
              break;

            case TileAlignment.BottomLeft:
              tileImageAbsolutePosition = vec(
                tileAbsolutePosition.x,
                tileAbsolutePosition.y + this.options.tileSize - tileImage.height
              );
              break;

            case TileAlignment.Bottom:
              tileImageAbsolutePosition = vec(
                (
                  tileAbsolutePosition.x + this.options.tileSize / 2
                ) - tileImage.width / 2,
                tileAbsolutePosition.y + this.options.tileSize - tileImage.height
              );
              break;

            case TileAlignment.BottomRight:
              tileImageAbsolutePosition = vec(
                tileAbsolutePosition.x + this.options.tileSize - tileImage.width,
                tileAbsolutePosition.y + this.options.tileSize - tileImage.height
              );
              break;
          }

          chunkContext.drawImage(
            tileImage,
            tileImageAbsolutePosition.x,
            tileImageAbsolutePosition.y
          );

          if (layer.clip) {
            chunkContext.restore();
          }

          layer.postRenderTile?.(
            chunkCanvas,
            chunkContext,
            this,
            layer,
            chunkPosition,
            tilePosition
          );
        }
      }

      chunkContext.restore();
    }

    this.options.postGenerateChunk?.(
      chunkCanvas,
      chunkContext,
      this,
      {
        topLeft: topLeftTile,
        bottomRight: bottomRightTile,
      },
      chunkPosition
    );

    return chunk;
  }

  private drawLine(
    context: CanvasRenderingContext2D,
    start: vec,
    end: vec,
    colour: string,
    lineWidth: number
  ) {
    context.save();

    context.lineWidth = lineWidth;
    context.strokeStyle = colour;

    context.beginPath();
    context.moveTo(start.x, start.y);
    context.lineTo(end.x, end.y);
    context.stroke();

    context.restore();
  }

  private drawCross(
    context: CanvasRenderingContext2D,
    position: vec,
    colour: string,
    lineWidth: number,
    size: number
  ) {
    context.save();

    context.lineWidth = lineWidth;

    const halfSize = Math.ceil(size / 2);
    context.strokeStyle = colour;
    context.beginPath();
    context.moveTo(position.x - halfSize, position.y - halfSize);
    context.lineTo(position.x + halfSize, position.y + halfSize);
    context.moveTo(position.x - halfSize, position.y + halfSize);
    context.lineTo(position.x + halfSize, position.y - halfSize);
    context.stroke();

    context.restore();
  }
}

/**
 * Content Manager Processor wrapper which converts TileMapOptionsData into
 * TileMapOptions
 *
 * @see https://www.npmjs.com/package/@basementuniverse/content-manager
 */
export async function tileMapOptionsContentProcessor(
  content: Record<string, {
    name: string;
    type: string;
    content: any;
    status: number;
  }>,
  data: {
    name: string;
    type: string;
    content: TileMapOptionsData;
    status: number;
  }
): Promise<void> {
  //
}
