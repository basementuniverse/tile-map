import { vec2 } from '@basementuniverse/vec';

export type Rectangle = {
  position: vec2;
  size: vec2;
};

export function bitmapToRectangles(bitmap: boolean[][]): Rectangle[] {
  const rectangles: Rectangle[] = [];

  // Step 1 - create 1-unit tall rectangles for each row
  for (const [y, row] of bitmap.entries()) {
    let currentRectangle: Rectangle | null = null;

    for (let x = 0; x < row.length; x++) {
      if (row[x]) {
        if (!currentRectangle) {
          currentRectangle = {
            position: vec2(x, y),
            size: vec2(1, 1),
          };
        } else {
          currentRectangle.size.x++;
        }
      } else {
        if (currentRectangle) {
          rectangles.push(currentRectangle);
          currentRectangle = null;
        }
      }
    }
  }

  // Step 2 - extend each rectangle downwards if possible
  let pair: [Rectangle, Rectangle] | null;
  while ((pair = findRectangleToExtend(rectangles))) {
    const [a, b] = pair;

    rectangles.splice(indexOf(b, rectangles), 1, ...chopRectangle(b, a));

    a.size.y += b.size.y;
  }

  return rectangles;
}

/**
 * Get the index of rectangle a in a list of rectangles
 */
function indexOf(a: Rectangle, rectangles: Rectangle[]) {
  return rectangles.findIndex(
    b => vec2.eq(a.position, b.position) && vec2.eq(a.size, b.size)
  );
}

/**
 * Find a pair of rectangles where the first one can be extended into the
 * second one
 *
 * If no such pair exists, return null
 */
function findRectangleToExtend(
  rectangles: Rectangle[]
): [Rectangle, Rectangle] | null {
  for (const a of rectangles) {
    const b = findRectangleToExtendInto(a, rectangles);
    if (b) {
      return [a, b];
    }
  }

  return null;
}

/**
 * Find a rectangle which rectangle a can extend into, or null if none can be
 * found
 *
 * A rectangle can extend into another one if the other one is exactly below
 * and their x-axis projections overlap
 */
function findRectangleToExtendInto(
  a: Rectangle,
  rectangles: Rectangle[]
): Rectangle | null {
  return (
    rectangles.find(
      other =>
        // The other rectangle is exactly below the current one
        other.position.y === a.position.y + a.size.y &&
        // The other rectangle starts before (or at) the start of the current one
        other.position.x <= a.position.x &&
        // The other rectangle ends after (or at) the end of the current one
        other.position.x + other.size.x >= a.position.x + a.size.x
    ) ?? null
  );
}

/**
 * Subtract rectangle b from rectangle a, ignoring height (i.e. only in the
 * x-axis) and return 0, 1 or 2 resulting rectangles
 */
function chopRectangle(a: Rectangle, b: Rectangle): Rectangle[] {
  const result: Rectangle[] = [];
  if (b.position.x > a.position.x) {
    result.push({
      position: vec2(a.position.x, a.position.y),
      size: vec2(b.position.x - a.position.x, a.size.y),
    });
  }
  if (b.position.x + b.size.x < a.position.x + a.size.x) {
    result.push({
      position: vec2(b.position.x + b.size.x, a.position.y),
      size: vec2(a.position.x + a.size.x - (b.position.x + b.size.x), a.size.y),
    });
  }

  return result;
}
