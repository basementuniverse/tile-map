(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define([], factory);
	else {
		var a = factory();
		for(var i in a) (typeof exports === 'object' ? exports : root)[i] = a[i];
	}
})(self, () => {
return /******/ (() => { // webpackBootstrap
/******/ 	var __webpack_modules__ = ({

/***/ "./node_modules/@basementuniverse/utils/utils.js":
/*!*******************************************************!*\
  !*** ./node_modules/@basementuniverse/utils/utils.js ***!
  \*******************************************************/
/***/ ((module) => {

/**
 * @overview A library of useful functions
 * @author Gordon Larrigan
 */

/**
 * Check if two numbers are approximately equal
 * @param {number} a Number a
 * @param {number} b Number b
 * @param {number} [p=Number.EPSILON] The precision value
 * @return {boolean} True if numbers a and b are approximately equal
 */
const floatEquals = (a, b, p = Number.EPSILON) => Math.abs(a - b) < p;

/**
 * Clamp a number between min and max
 * @param {number} a The number to clamp
 * @param {number} [min=0] The minimum value
 * @param {number} [max=1] The maximum value
 * @return {number} A clamped number
 */
const clamp = (a, min = 0, max = 1) => a < min ? min : (a > max ? max : a);

/**
 * Get the fractional part of a number
 * @param {number} a The number from which to get the fractional part
 * @return {number} The fractional part of the number
 */
const frac = a => a >= 0 ? a - Math.floor(a) : a - Math.ceil(a);

/**
 * Round n to d decimal places
 * @param {number} n The number to round
 * @param {number} [d=0] The number of decimal places to round to
 * @return {number} A rounded number
 */
const round = (n, d = 0) => {
  const p = Math.pow(10, d);
  return Math.round(n * p + Number.EPSILON) / p;
}

/**
 * Do a linear interpolation between a and b
 * @param {number} a The minimum number
 * @param {number} b The maximum number
 * @param {number} i The interpolation value, should be in the interval [0, 1]
 * @return {number} An interpolated value in the interval [a, b]
 */
const lerp = (a, b, i) => a + (b - a) * i;

/**
 * Get the position of i between a and b
 * @param {number} a The minimum number
 * @param {number} b The maximum number
 * @param {number} i The interpolated value in the interval [a, b]
 * @return {number} The position of i between a and b
 */
const unlerp = (a, b, i) => (i - a) / (b - a);

/**
 * Do a bilinear interpolation
 * @param {number} c00 Top-left value
 * @param {number} c10 Top-right value
 * @param {number} c01 Bottom-left value
 * @param {number} c11 Bottom-right value
 * @param {number} ix Interpolation value along x
 * @param {number} iy Interpolation value along y
 * @return {number} A bilinear interpolated value
 */
const blerp = (c00, c10, c01, c11, ix, iy) => lerp(lerp(c00, c10, ix), lerp(c01, c11, ix), iy);

/**
 * Re-map a number i from range a1...a2 to b1...b2
 * @param {number} i The number to re-map
 * @param {number} a1
 * @param {number} a2
 * @param {number} b1
 * @param {number} b2
 * @return {number}
 */
const remap = (i, a1, a2, b1, b2) => b1 + (i - a1) * (b2 - b1) / (a2 - a1);

/**
 * Do a smooth interpolation between a and b
 * @param {number} a The minimum number
 * @param {number} b The maximum number
 * @param {number} i The interpolation value
 * @return {number} An interpolated value in the interval [a, b]
 */
const smoothstep = (a, b, i) => lerp(a, b, 3 * Math.pow(i, 2) - 2 * Math.pow(i, 3));

/**
 * Get an angle in radians
 * @param {number} degrees The angle in degrees
 * @return {number} The angle in radians
 */
const radians = degrees => (Math.PI / 180) * degrees;

/**
 * Get an angle in degrees
 * @param {number} radians The angle in radians
 * @return {number} The angle in degrees
 */
const degrees = radians => (180 / Math.PI) * radians;

/**
 * Get a random float in the interval [min, max)
 * @param {number} min Inclusive min
 * @param {number} max Exclusive max
 * @return {number} A random float in the interval [min, max)
 */
const randomBetween = (min, max) => Math.random() * (max - min) + min;

/**
 * Get a random integer in the interval [min, max]
 * @param {number} min Inclusive min
 * @param {number} max Inclusive max
 * @return {number} A random integer in the interval [min, max]
 */
const randomIntBetween = (min, max) => Math.floor(Math.random() * (max - min + 1)) + min;

/**
 * Get a normally-distributed random number
 * @param {number} [mu=0.5] The mean value
 * @param {number} [sigma=0.5] The standard deviation
 * @param {number} [samples=2] The number of samples
 * @return {number} A normally-distributed random number
 */
const cltRandom = (mu = 0.5, sigma = 0.5, samples = 2) => {
  let total = 0;
  for (let i = samples; i--;) {
    total += Math.random();
  }
  return mu + (total - samples / 2) / (samples / 2) * sigma;
};

/**
 * Get a normally-distributed random integer in the interval [min, max]
 * @param {number} min Inclusive min
 * @param {number} max Inclusive max
 * @return {number} A normally-distributed random integer
 */
const cltRandomInt = (min, max) => Math.floor(min + cltRandom(0.5, 0.5, 2) * (max + 1 - min));

/**
 * Return a weighted random integer
 * @param {Array<number>} w An array of weights
 * @return {number} An index from w
 */
const weightedRandom = w => {
  let total = w.reduce((a, i) => a + i, 0), n = 0;
  const r = Math.random() * total;
  while (total > r) {
    total -= w[n++];
  }
  return n - 1;
};

/**
 * An interpolation function
 * @callback InterpolationFunction
 * @param {number} a The minimum number
 * @param {number} b The maximum number
 * @param {number} i The interpolation value, should be in the interval [0, 1]
 * @return {number} The interpolated value in the interval [a, b]
 */

/**
 * Return an interpolated value from an array
 * @param {Array<number>} a An array of values interpolate
 * @param {number} i A number in the interval [0, 1]
 * @param {InterpolationFunction} [f=Math.lerp] The interpolation function to use
 * @return {number} An interpolated value in the interval [min(a), max(a)]
 */
const lerpArray = (a, i, f = lerp) => {
  const s = i * (a.length - 1);
  const p = clamp(Math.trunc(s), 0, a.length - 1);
  return f(a[p] || 0, a[p + 1] || 0, frac(s));
};

/**
 * Get the dot product of two vectors
 * @param {Array<number>} a Vector a
 * @param {Array<number>} b Vector b
 * @return {number} a ∙ b
 */
const dot = (a, b) => a.reduce((n, v, i) => n + v * b[i], 0);

/**
 * Get the factorial of a number
 * @param {number} a
 * @return {number} a!
 */
const factorial = a => {
  let result = 1;
  for (let i = 2; i <= a; i++) {
    result *= i;
  }
  return result;
};

/**
 * Get the number of permutations of r elements from a set of n elements
 * @param {number} n
 * @param {number} r
 * @return {number} nPr
 */
const npr = (n, r) => factorial(n) / factorial(n - r);

/**
 * Get the number of combinations of r elements from a set of n elements
 * @param {number} n
 * @param {number} r
 * @return {number} nCr
 */
const ncr = (n, r) => factorial(n) / (factorial(r) * factorial(n - r));

/**
 * Generate all combinations of r elements from an array
 *
 * @example
 * ```js
 * combinations([1, 2, 3], 2);
 * ```
 *
 * Output:
 * ```json
 * [
 *   [1, 2],
 *   [1, 3],
 *   [2, 3]
 * ]
 * ```
 * @param {Array<*>} a
 * @param {number} r The number of elements to choose in each combination
 * @return {Array<Array<*>>} An array of combination arrays
 */
const combinations = (a, r) => {
  if (r === 1) {
    return a.map(item => [item]);
  }

  return a.reduce(
    (acc, item, i) => [
      ...acc,
      ...combinations(a.slice(i + 1), r - 1).map(c => [item, ...c]),
    ],
    []
  );
};

/**
 * Get a cartesian product of arrays
 *
 * @example
 * ```js
 * cartesian([1, 2, 3], ['a', 'b']);
 * ```
 *
 * Output:
 * ```json
 * [
 *   [1, "a"],
 *   [1, "b"],
 *   [2, "a"],
 *   [2, "b"],
 *   [3, "a"],
 *   [3, "b"]
 * ]
 * ```
 */
const cartesian = (...arr) =>
  arr.reduce(
    (a, b) => a.flatMap(c => b.map(d => [...c, d])),
    [[]]
  );

/**
 * A function for generating array values
 * @callback TimesFunction
 * @param {number} i The array index
 * @return {*} The array value
 */

/**
 * Return a new array with length n by calling function f(i) on each element
 * @param {TimesFunction} f
 * @param {number} n The size of the array
 * @return {Array<*>}
 */
const times = (f, n) => Array(n).fill(0).map((_, i) => f(i));

/**
 * Return an array containing numbers 0->(n - 1)
 * @param {number} n The size of the array
 * @return {Array<number>} An array of integers 0->(n - 1)
 */
const range = n => times(i => i, n);

/**
 * Zip 2 arrays together, i.e. ([1, 2, 3], [a, b, c]) => [[1, a], [2, b], [3, c]]
 * @param {Array<*>} a
 * @param {Array<*>} b
 * @return {Array<Array<*>>}
 */
const zip = (a, b) => a.map((k, i) => [k, b[i]]);

/**
 * Return array[i] with positive and negative wrapping
 * @param {Array<*>} a
 * @param {number} i The positively/negatively wrapped array index
 * @return {*} An element from the array
 */
const at = (a, i) => a[i < 0 ? a.length - (Math.abs(i + 1) % a.length) - 1 : i % a.length];

/**
 * Return the last element of an array without removing it
 * @param {Array<*>} a
 * @return {*} The last element from the array
 */
const peek = (a) => {
  if (!a.length) {
    return undefined;
  }

  return a[a.length - 1];
};

/**
 * Chop an array into chunks of size n
 * @param {Array<*>} a
 * @param {number} n The chunk size
 * @return {Array<Array<*>>} An array of array chunks
 */
const chunk = (a, n) => times(i => a.slice(i * n, i * n + n), Math.ceil(a.length / n));

/**
 * Randomly shuffle a shallow copy of an array
 * @param {Array<*>} a
 * @return {Array<*>} The shuffled array
 */
const shuffle = a => a.slice().sort(() => Math.random() - 0.5);

/**
 * Flatten an object
 * @param {object} o
 * @param {string} concatenator The string to use for concatenating keys
 * @return {object} A flattened object
 */
const flat = (o, concatenator = '.') => {
  return Object.keys(o).reduce((acc, key) => {
    if (o[key] instanceof Date) {
      return {
        ...acc,
        [key]: o[key].toISOString(),
      };
    }

    if (typeof o[key] !== 'object' || !o[key]) {
      return {
        ...acc,
        [key]: o[key],
      };
    }
    const flattened = flat(o[key], concatenator);

    return {
      ...acc,
      ...Object.keys(flattened).reduce(
        (childAcc, childKey) => ({
          ...childAcc,
          [`${key}${concatenator}${childKey}`]: flattened[childKey],
        }),
        {}
      ),
    };
  }, {});
};

/**
 * Unflatten an object
 * @param {object} o
 * @param {string} concatenator The string to check for in concatenated keys
 * @return {object} An un-flattened object
 */
const unflat = (o, concatenator = '.') => {
  let result = {}, temp, substrings, property, i;

  for (property in o) {
    substrings = property.split(concatenator);
    temp = result;
    for (i = 0; i < substrings.length - 1; i++) {
      if (!(substrings[i] in temp)) {
        if (isFinite(substrings[i + 1])) {
          temp[substrings[i]] = [];
        } else {
          temp[substrings[i]] = {};
        }
      }
      temp = temp[substrings[i]];
    }
    temp[substrings[substrings.length - 1]] = o[property];
  }

  return result;
};

/**
 * A split predicate
 * @callback SplitPredicate
 * @param {any} value The current value
 * @return {boolean} True if the array should split at this index
 */

/**
 * Split an array into sub-arrays based on a predicate
 * @param {Array<*>} array
 * @param {SplitPredicate} predicate
 * @return {Array<Array<*>>} An array of arrays
 */
const split = (array, predicate) => {
  const result = [];
  let current = [];
  for (const value of array) {
    if (predicate(value)) {
      if (current.length) {
        result.push(current);
      }
      current = [value];
    } else {
      current.push(value);
    }
  }
  result.push(current);

  return result;
};

/**
 * Pluck keys from an object
 * @param {object} o
 * @param {...string} keys The keys to pluck from the object
 * @return {object} An object containing the plucked keys
 */
const pluck = (o, ...keys) => {
  return keys.reduce(
    (result, key) => Object.assign(result, { [key]: o[key] }),
    {}
  );
};

/**
 * Exclude keys from an object
 * @param {object} o
 * @param {...string} keys The keys to exclude from the object
 * @return {object} An object containing all keys except excluded keys
 */
const exclude = (o, ...keys) => {
  return Object.fromEntries(
    Object.entries(o).filter(([key]) => !keys.includes(key))
  );
};

if (true) {
  module.exports = {
    floatEquals,
    clamp,
    frac,
    round,
    lerp,
    unlerp,
    blerp,
    remap,
    smoothstep,
    radians,
    degrees,
    randomBetween,
    randomIntBetween,
    cltRandom,
    cltRandomInt,
    weightedRandom,
    lerpArray,
    dot,
    factorial,
    npr,
    ncr,
    combinations,
    cartesian,
    times,
    range,
    zip,
    at,
    peek,
    chunk,
    shuffle,
    flat,
    unflat,
    split,
    pluck,
    exclude,
  };
}


/***/ }),

/***/ "./node_modules/@basementuniverse/vec/vec.js":
/*!***************************************************!*\
  !*** ./node_modules/@basementuniverse/vec/vec.js ***!
  \***************************************************/
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {

const { times, chunk, dot } = __webpack_require__(/*! @basementuniverse/utils */ "./node_modules/@basementuniverse/utils/utils.js");

/**
 * @overview A small vector and matrix library
 * @author Gordon Larrigan
 */

/**
 * A 2d vector
 * @typedef {Object} vec
 * @property {number} x The x component of the vector
 * @property {number} y The y component of the vector
 */

/**
 * Create a new vector
 * @param {number|vec} [x] The x component of the vector, or a vector to copy
 * @param {number} [y] The y component of the vector
 * @return {vec} A new vector
 * @example <caption>Various ways to initialise a vector</caption>
 * let a = vec(3, 2);  // (3, 2)
 * let b = vec(4);     // (4, 4)
 * let c = vec(a);     // (3, 2)
 * let d = vec();      // (0, 0)
 */
const vec = (x, y) => (!x && !y ?
  { x: 0, y: 0 } : (typeof x === 'object' ?
    { x: x.x || 0, y: x.y || 0 } : (y === null || y === undefined ?
      { x: x, y: x } : { x: x, y: y })
  )
);

/**
 * Get the components of a vector as an array
 * @param {vec} a The vector to get components from
 * @return {Array<number>} The vector components as an array
 */
vec.components = a => [a.x, a.y];

/**
 * Return a unit vector (1, 0)
 * @return {vec} A unit vector (1, 0)
 */
vec.ux = () => vec(1, 0);

/**
 * Return a unit vector (0, 1)
 * @return {vec} A unit vector (0, 1)
 */
vec.uy = () => vec(0, 1);

/**
 * Add vectors
 * @param {vec} a Vector a
 * @param {vec} b Vector b
 * @return {vec} a + b
 */
vec.add = (a, b) => ({ x: a.x + b.x, y: a.y + b.y });

/**
 * Scale a vector
 * @param {vec} a Vector a
 * @param {number} b Scalar b
 * @return {vec} a * b
 */
vec.mul = (a, b) => ({ x: a.x * b, y: a.y * b });

/**
 * Subtract vectors
 * @param {vec} a Vector a
 * @param {vec} b Vector b
 * @return {vec} a - b
 */
vec.sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });

/**
 * Get the length of a vector
 * @param {vec} a Vector a
 * @return {number} |a|
 */
vec.len = a => Math.sqrt(a.x * a.x + a.y * a.y);

/**
 * Get the length of a vector using taxicab geometry
 * @param {vec} a Vector a
 * @return {number} |a|
 */
vec.manhattan = a => Math.abs(a.x) + Math.abs(a.y);

/**
 * Normalise a vector
 * @param {vec} a The vector to normalise
 * @return {vec} ^a
 */
vec.nor = a => {
  let len = vec.len(a);
  return len ? { x: a.x / len, y: a.y / len } : vec();
};

/**
 * Get a dot product of vectors
 * @param {vec} a Vector a
 * @param {vec} b Vector b
 * @return {number} a ∙ b
 */
vec.dot = (a, b) => a.x * b.x + a.y * b.y;

/**
 * Rotate a vector by r radians
 * @param {vec} a The vector to rotate
 * @param {number} r The angle to rotate by, measured in radians
 * @return {vec} A rotated vector
 */
vec.rot = (a, r) => {
  let s = Math.sin(r),
    c = Math.cos(r);
  return { x: c * a.x - s * a.y, y: s * a.x + c * a.y };
}

/**
 * Check if two vectors are equal
 * @param {vec} a Vector a
 * @param {vec} b Vector b
 * @return {boolean} True if vectors a and b are equal, false otherwise
 */
vec.eq = (a, b) => a.x === b.x && a.y === b.y;

/**
 * Get the angle of a vector
 * @param {vec} a Vector a
 * @return {number} The angle of vector a in radians
 */
vec.rad = a => Math.atan2(a.y, a.x);

/**
 * Copy a vector
 * @param {vec} a The vector to copy
 * @return {vec} A copy of vector a
 */
vec.cpy = a => vec(a);

/**
 * A function to call on each component of a vector
 * @callback vectorMapCallback
 * @param {number} value The component value
 * @param {'x' | 'y'} label The component label (x or y)
 * @return {number} The mapped component
 */

/**
 * Call a function on each component of a vector and build a new vector from the results
 * @param {vec} a Vector a
 * @param {vectorMapCallback} f The function to call on each component of the vector
 * @return {vec} Vector a mapped through f
 */
vec.map = (a, f) => ({ x: f(a.x, 'x'), y: f(a.y, 'y') });

/**
 * Convert a vector into a string
 * @param {vec} a The vector to convert
 * @param {string} [s=', '] The separator string
 * @return {string} A string representation of the vector
 */
vec.str = (a, s = ', ') => `${a.x}${s}${a.y}`;

/**
 * A matrix
 * @typedef {Object} mat
 * @property {number} m The number of rows in the matrix
 * @property {number} n The number of columns in the matrix
 * @property {Array<number>} entries The matrix values
 */

/**
 * Create a new matrix
 * @param {number} [m=4] The number of rows
 * @param {number} [n=4] The number of columns
 * @param {Array<number>} [entries=[]] Matrix values in reading order
 * @return {mat} A new matrix
 */
const mat = (m = 4, n = 4, entries = []) => ({
  m, n,
  entries: entries.concat(Array(m * n).fill(0)).slice(0, m * n)
});

/**
 * Get an identity matrix of size n
 * @param {number} n The size of the matrix
 * @return {mat} An identity matrix
 */
mat.identity = n => mat(n, n, Array(n * n).fill(0).map((v, i) => +(Math.floor(i / n) === i % n)));

/**
 * Get an entry from a matrix
 * @param {mat} a Matrix a
 * @param {number} i The row offset
 * @param {number} j The column offset
 * @return {number} The value at position (i, j) in matrix a
 */
mat.get = (a, i, j) => a.entries[(j - 1) + (i - 1) * a.n];

/**
 * Set an entry of a matrix
 * @param {mat} a Matrix a
 * @param {number} i The row offset
 * @param {number} j The column offset
 * @param {number} v The value to set in matrix a
 */
mat.set = (a, i, j, v) => { a.entries[(j - 1) + (i - 1) * a.n] = v; };

/**
 * Get a row from a matrix as an array
 * @param {mat} a Matrix a
 * @param {number} m The row offset
 * @return {Array<number>} Row m from matrix a
 */
mat.row = (a, m) => {
  const s = (m - 1) * a.n;
  return a.entries.slice(s, s + a.n);
};

/**
 * Get a column from a matrix as an array
 * @param {mat} a Matrix a
 * @param {number} n The column offset
 * @return {Array<number>} Column n from matrix a
 */
mat.col = (a, n) => times(i => mat.get(a, (i + 1), n), a.m);

/**
 * Add matrices
 * @param {mat} a Matrix a
 * @param {mat} b Matrix b
 * @return {mat} a + b
 */
mat.add = (a, b) => a.m === b.m && a.n === b.n && mat.map(a, (v, i) => v + b.entries[i]);

/**
 * Subtract matrices
 * @param {mat} a Matrix a
 * @param {mat} b Matrix b
 * @return {mat} a - b
 */
mat.sub = (a, b) => a.m === b.m && a.n === b.n && mat.map(a, (v, i) => v - b.entries[i]);

/**
 * Multiply matrices
 * @param {mat} a Matrix a
 * @param {mat} b Matrix b
 * @return {mat|boolean} ab or false if the matrices cannot be multiplied
 */
mat.mul = (a, b) => {
  if (a.n !== b.m) { return false; }
  const result = mat(a.m, b.n);
  for (let i = 1; i <= a.m; i++) {
    for (let j = 1; j <= b.n; j++) {
      mat.set(result, i, j, dot(mat.row(a, i), mat.col(b, j)));
    }
  }
  return result;
};

/**
 * Scale a matrix
 * @param {mat} a Matrix a
 * @param {number} b Scalar b
 * @return {mat} a * b
 */
mat.scale = (a, b) => mat.map(a, v => v * b);

/**
 * Transpose a matrix
 * @param {mat} a The matrix to transpose
 * @return {mat} A transposed matrix
 */
mat.trans = a => mat(a.n, a.m, times(i => mat.col(a, (i + 1)), a.n).flat());

/**
 * Get the minor of a matrix
 * @param {mat} a Matrix a
 * @param {number} i The row offset
 * @param {number} j The column offset
 * @return {mat|boolean} The (i, j) minor of matrix a or false if the matrix is not square
 */
mat.minor = (a, i, j) => {
  if (a.m !== a.n) { return false; }
  const entries = [];
  for (let ii = 1; ii <= a.m; ii++) {
    if (ii === i) { continue; }
    for (let jj = 1; jj <= a.n; jj++) {
      if (jj === j) { continue; }
      entries.push(mat.get(a, ii, jj));
    }
  }
  return mat(a.m - 1, a.n - 1, entries);
};

/**
 * Get the determinant of a matrix
 * @param {mat} a Matrix a
 * @return {number|boolean} |a| or false if the matrix is not square
 */
mat.det = a => {
  if (a.m !== a.n) { return false; }
  if (a.m === 1) {
    return a.entries[0];
  }
  if (a.m === 2) {
    return a.entries[0] * a.entries[3] - a.entries[1] * a.entries[2];
  }
  let total = 0, sign = 1;
  for (let j = 1; j <= a.n; j++) {
    total += sign * a.entries[j - 1] * mat.det(mat.minor(a, 1, j));
    sign *= -1;
  }
  return total;
};

/**
 * Normalise a matrix
 * @param {mat} a The matrix to normalise
 * @return {mat|boolean} ^a or false if the matrix is not square
 */
mat.nor = a => {
  if (a.m !== a.n) { return false; }
  const d = mat.det(a);
  return mat.map(a, i => i * d);
};

/**
 * Get the adjugate of a matrix
 * @param {mat} a The matrix from which to get the adjugate
 * @return {mat} The adjugate of a
 */
mat.adj = a => {
  const minors = mat(a.m, a.n);
  for (let i = 1; i <= a.m; i++) {
    for (let j = 1; j <= a.n; j++) {
      mat.set(minors, i, j, mat.det(mat.minor(a, i, j)));
    }
  }
  const cofactors = mat.map(minors, (v, i) => v * (i % 2 ? -1 : 1));
  return mat.trans(cofactors);
};

/**
 * Get the inverse of a matrix
 * @param {mat} a The matrix to invert
 * @return {mat|boolean} a^-1 or false if the matrix has no inverse
 */
mat.inv = a => {
  if (a.m !== a.n) { return false; }
  const d = mat.det(a);
  if (d === 0) { return false; }
  return mat.scale(mat.adj(a), 1 / d);
};

/**
 * Check if two matrices are equal
 * @param {mat} a Matrix a
 * @param {mat} b Matrix b
 * @return {boolean} True if matrices a and b are identical, false otherwise
 */
mat.eq = (a, b) => a.m === b.m && a.n === b.n && mat.str(a) === mat.str(b);

/**
 * Copy a matrix
 * @param {mat} a The matrix to copy
 * @return {mat} A copy of matrix a
 */
mat.cpy = a => mat(a.m, a.n, [...a.entries]);

/**
 * A function to call on each entry of a matrix
 * @callback matrixMapCallback
 * @param {number} value The entry value
 * @param {number} index The entry index
 * @param {Array<number>} entries The array of matrix entries
 * @return {number} The mapped entry
 */

/**
 * Call a function on each entry of a matrix and build a new matrix from the results
 * @param {mat} a Matrix a
 * @param {matrixMapCallback} f The function to call on each entry of the matrix
 * @return {mat} Matrix a mapped through f
 */
mat.map = (a, f) => mat(a.m, a.n, a.entries.map(f));

/**
 * Convert a matrix into a string
 * @param {mat} a The matrix to convert
 * @param {string} [ms=', '] The separator string for columns
 * @param {string} [ns='\n'] The separator string for rows
 * @return {string} A string representation of the matrix
 */
mat.str = (a, ms = ', ', ns = '\n') => chunk(a.entries, a.n).map(r => r.join(ms)).join(ns);

if (true) {
  module.exports = { vec, mat };
}


/***/ }),

/***/ "./node_modules/fast-rle/decode.js":
/*!*****************************************!*\
  !*** ./node_modules/fast-rle/decode.js ***!
  \*****************************************/
/***/ ((module) => {

const decode = nums => {
  const decoded = [];
  for (let i = 0; i < nums.length; i += 2) {
    const run_length = nums[i];
    const value = nums[i + 1];
    for (let ii = 0; ii < run_length; ii++) {
      decoded.push(value);
    }
  }
  return decoded;
};

module.exports = decode;


/***/ }),

/***/ "./node_modules/lru_map/dist/lru.js":
/*!******************************************!*\
  !*** ./node_modules/lru_map/dist/lru.js ***!
  \******************************************/
/***/ (function(__unused_webpack_module, exports) {

!function(g,c){ true?c(exports):0}(this,function(g){const c=Symbol("newer"),e=Symbol("older");class n{constructor(a,b){typeof a!=="number"&&(b=a,a=0),this.size=0,this.limit=a,this.oldest=this.newest=void 0,this._keymap=new Map(),b&&(this.assign(b),a<1&&(this.limit=this.size))}_markEntryAsUsed(a){if(a===this.newest)return;a[c]&&(a===this.oldest&&(this.oldest=a[c]),a[c][e]=a[e]),a[e]&&(a[e][c]=a[c]),a[c]=void 0,a[e]=this.newest,this.newest&&(this.newest[c]=a),this.newest=a}assign(a){let b,d=this.limit||Number.MAX_VALUE;this._keymap.clear();let m=a[Symbol.iterator]();for(let h=m.next();!h.done;h=m.next()){let f=new l(h.value[0],h.value[1]);this._keymap.set(f.key,f),b?(b[c]=f,f[e]=b):this.oldest=f,b=f;if(d--==0)throw new Error("overflow")}this.newest=b,this.size=this._keymap.size}get(a){var b=this._keymap.get(a);return b?(this._markEntryAsUsed(b),b.value):void 0}set(a,b){var d=this._keymap.get(a);return d?(d.value=b,this._markEntryAsUsed(d),this):(this._keymap.set(a,d=new l(a,b)),this.newest?(this.newest[c]=d,d[e]=this.newest):this.oldest=d,this.newest=d,++this.size,this.size>this.limit&&this.shift(),this)}shift(){var a=this.oldest;if(a)return this.oldest[c]?(this.oldest=this.oldest[c],this.oldest[e]=void 0):(this.oldest=void 0,this.newest=void 0),a[c]=a[e]=void 0,this._keymap.delete(a.key),--this.size,[a.key,a.value]}find(a){let b=this._keymap.get(a);return b?b.value:void 0}has(a){return this._keymap.has(a)}delete(a){var b=this._keymap.get(a);return b?(this._keymap.delete(b.key),b[c]&&b[e]?(b[e][c]=b[c],b[c][e]=b[e]):b[c]?(b[c][e]=void 0,this.oldest=b[c]):b[e]?(b[e][c]=void 0,this.newest=b[e]):this.oldest=this.newest=void 0,this.size--,b.value):void 0}clear(){this.oldest=this.newest=void 0,this.size=0,this._keymap.clear()}keys(){return new j(this.oldest)}values(){return new k(this.oldest)}entries(){return this}[Symbol.iterator](){return new i(this.oldest)}forEach(a,b){typeof b!=="object"&&(b=this);let d=this.oldest;for(;d;)a.call(b,d.value,d.key,this),d=d[c]}toJSON(){for(var a=new Array(this.size),b=0,d=this.oldest;d;)a[b++]={key:d.key,value:d.value},d=d[c];return a}toString(){for(var a="",b=this.oldest;b;)a+=String(b.key)+":"+b.value,b=b[c],b&&(a+=" < ");return a}}g.LRUMap=n;function l(a,b){this.key=a,this.value=b,this[c]=void 0,this[e]=void 0}function i(a){this.entry=a}i.prototype[Symbol.iterator]=function(){return this},i.prototype.next=function(){let a=this.entry;return a?(this.entry=a[c],{done:!1,value:[a.key,a.value]}):{done:!0,value:void 0}};function j(a){this.entry=a}j.prototype[Symbol.iterator]=function(){return this},j.prototype.next=function(){let a=this.entry;return a?(this.entry=a[c],{done:!1,value:a.key}):{done:!0,value:void 0}};function k(a){this.entry=a}k.prototype[Symbol.iterator]=function(){return this},k.prototype.next=function(){let a=this.entry;return a?(this.entry=a[c],{done:!1,value:a.value}):{done:!0,value:void 0}}});
//# sourceMappingURL=lru.js.map


/***/ }),

/***/ "./bitmap-decompose.ts":
/*!*****************************!*\
  !*** ./bitmap-decompose.ts ***!
  \*****************************/
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";

Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.bitmapToRectangles = void 0;
const vec_1 = __webpack_require__(/*! @basementuniverse/vec */ "./node_modules/@basementuniverse/vec/vec.js");
function bitmapToRectangles(bitmap) {
    const rectangles = [];
    // Step 1 - create 1-unit tall rectangles for each row
    for (const [y, row] of bitmap.entries()) {
        let currentRectangle = null;
        for (let x = 0; x < row.length; x++) {
            if (row[x]) {
                if (!currentRectangle) {
                    currentRectangle = {
                        position: (0, vec_1.vec)(x, y),
                        size: (0, vec_1.vec)(1, 1),
                    };
                }
                else {
                    currentRectangle.size.x++;
                }
            }
            else {
                if (currentRectangle) {
                    rectangles.push(currentRectangle);
                    currentRectangle = null;
                }
            }
        }
    }
    // Step 2 - extend each rectangle downwards if possible
    let pair;
    while (pair = findRectangleToExtend(rectangles)) {
        const [a, b] = pair;
        rectangles.splice(indexOf(b, rectangles), 1, ...chopRectangle(b, a));
        a.size.y += b.size.y;
    }
    return rectangles;
}
exports.bitmapToRectangles = bitmapToRectangles;
/**
 * Get the index of rectangle a in a list of rectangles
 */
function indexOf(a, rectangles) {
    return rectangles.findIndex(b => vec_1.vec.eq(a.position, b.position) && vec_1.vec.eq(a.size, b.size));
}
/**
 * Find a pair of rectangles where the first one can be extended into the
 * second one
 *
 * If no such pair exists, return null
 */
function findRectangleToExtend(rectangles) {
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
function findRectangleToExtendInto(a, rectangles) {
    var _a;
    return (_a = rectangles.find(other => (
    // The other rectangle is exactly below the current one
    other.position.y === a.position.y + a.size.y &&
        // The other rectangle starts before (or at) the start of the current one
        other.position.x <= a.position.x &&
        // The other rectangle ends after (or at) the end of the current one
        other.position.x + other.size.x >= a.position.x + a.size.x))) !== null && _a !== void 0 ? _a : null;
}
/**
 * Subtract rectangle b from rectangle a, ignoring height (i.e. only in the
 * x-axis) and return 0, 1 or 2 resulting rectangles
 */
function chopRectangle(a, b) {
    const result = [];
    if (b.position.x > a.position.x) {
        result.push({
            position: (0, vec_1.vec)(a.position.x, a.position.y),
            size: (0, vec_1.vec)(b.position.x - a.position.x, a.size.y),
        });
    }
    if (b.position.x + b.size.x < a.position.x + a.size.x) {
        result.push({
            position: (0, vec_1.vec)(b.position.x + b.size.x, a.position.y),
            size: (0, vec_1.vec)((a.position.x + a.size.x) - (b.position.x + b.size.x), a.size.y),
        });
    }
    return result;
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiYml0bWFwLWRlY29tcG9zZS5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzIjpbIi4uL2JpdG1hcC1kZWNvbXBvc2UudHMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7O0FBQUEsK0NBQTRDO0FBTzVDLFNBQWdCLGtCQUFrQixDQUFDLE1BQW1CO0lBQ3BELE1BQU0sVUFBVSxHQUFnQixFQUFFLENBQUM7SUFFbkMsc0RBQXNEO0lBQ3RELEtBQUssTUFBTSxDQUFDLENBQUMsRUFBRSxHQUFHLENBQUMsSUFBSSxNQUFNLENBQUMsT0FBTyxFQUFFLEVBQUU7UUFDdkMsSUFBSSxnQkFBZ0IsR0FBcUIsSUFBSSxDQUFDO1FBRTlDLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxHQUFHLENBQUMsTUFBTSxFQUFFLENBQUMsRUFBRSxFQUFFO1lBQ25DLElBQUksR0FBRyxDQUFDLENBQUMsQ0FBQyxFQUFFO2dCQUNWLElBQUksQ0FBQyxnQkFBZ0IsRUFBRTtvQkFDckIsZ0JBQWdCLEdBQUc7d0JBQ2pCLFFBQVEsRUFBRSxJQUFBLFNBQUcsRUFBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDO3dCQUNuQixJQUFJLEVBQUUsSUFBQSxTQUFHLEVBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQztxQkFDaEIsQ0FBQztpQkFDSDtxQkFBTTtvQkFDTCxnQkFBZ0IsQ0FBQyxJQUFJLENBQUMsQ0FBQyxFQUFFLENBQUM7aUJBQzNCO2FBQ0Y7aUJBQU07Z0JBQ0wsSUFBSSxnQkFBZ0IsRUFBRTtvQkFDcEIsVUFBVSxDQUFDLElBQUksQ0FBQyxnQkFBZ0IsQ0FBQyxDQUFDO29CQUNsQyxnQkFBZ0IsR0FBRyxJQUFJLENBQUM7aUJBQ3pCO2FBQ0Y7U0FDRjtLQUNGO0lBRUQsdURBQXVEO0lBQ3ZELElBQUksSUFBbUMsQ0FBQztJQUN4QyxPQUFPLElBQUksR0FBRyxxQkFBcUIsQ0FBQyxVQUFVLENBQUMsRUFBRTtRQUMvQyxNQUFNLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FB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

/***/ }),

/***/ "./index.ts":
/*!******************!*\
  !*** ./index.ts ***!
  \******************/
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {

"use strict";

var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.tileMapOptionsContentProcessor = exports.TileMap = exports.TileAlignment = void 0;
const lru_map_1 = __webpack_require__(/*! lru_map */ "./node_modules/lru_map/dist/lru.js");
const decode_1 = __importDefault(__webpack_require__(/*! fast-rle/decode */ "./node_modules/fast-rle/decode.js"));
const vec_1 = __webpack_require__(/*! @basementuniverse/vec */ "./node_modules/@basementuniverse/vec/vec.js");
const utils_1 = __webpack_require__(/*! @basementuniverse/utils */ "./node_modules/@basementuniverse/utils/utils.js");
const bitmap_decompose_1 = __webpack_require__(/*! ./bitmap-decompose */ "./bitmap-decompose.ts");
var TileAlignment;
(function (TileAlignment) {
    TileAlignment[TileAlignment["TopLeft"] = 0] = "TopLeft";
    TileAlignment[TileAlignment["Top"] = 1] = "Top";
    TileAlignment[TileAlignment["TopRight"] = 2] = "TopRight";
    TileAlignment[TileAlignment["Left"] = 3] = "Left";
    TileAlignment[TileAlignment["Center"] = 4] = "Center";
    TileAlignment[TileAlignment["Right"] = 5] = "Right";
    TileAlignment[TileAlignment["BottomLeft"] = 6] = "BottomLeft";
    TileAlignment[TileAlignment["Bottom"] = 7] = "Bottom";
    TileAlignment[TileAlignment["BottomRight"] = 8] = "BottomRight";
})(TileAlignment = exports.TileAlignment || (exports.TileAlignment = {}));
function pointInRectangle(point, topLeft, bottomRight) {
    return (point.x >= topLeft.x &&
        point.y >= topLeft.y &&
        point.x < bottomRight.x &&
        point.y < bottomRight.y);
}
class TileMap {
    constructor(options) {
        const actualOptions = Object.assign({}, TileMap.DEFAULT_OPTIONS, options !== null && options !== void 0 ? options : {});
        if (!actualOptions.debug || actualOptions.debug === true) {
            actualOptions.debug = {
                showOrigin: !!actualOptions.debug,
                showChunkBorders: !!actualOptions.debug,
                showChunkLabels: !!actualOptions.debug,
                showTileBorders: !!actualOptions.debug,
            };
        }
        this.options = actualOptions;
        this.chunkBuffer = new lru_map_1.LRUMap(this.options.chunkBufferMaxSize);
    }
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
    getLayerRectangles(layerName, fieldName, tileBounds) {
        var _a, _b, _c, _d, _e, _f, _g, _h, _j;
        const layer = this.options.layers.find((l) => l.name === layerName);
        if (!layer) {
            return [];
        }
        const topLeft = (_a = tileBounds === null || tileBounds === void 0 ? void 0 : tileBounds.topLeft) !== null && _a !== void 0 ? _a : (0, vec_1.vec)(0);
        const bottomRight = (_b = tileBounds === null || tileBounds === void 0 ? void 0 : tileBounds.bottomRight) !== null && _b !== void 0 ? _b : (0, vec_1.vec)(Math.max(...(_d = (_c = layer.data) === null || _c === void 0 ? void 0 : _c.map(row => row.length)) !== null && _d !== void 0 ? _d : [0]), (_f = (_e = layer.data) === null || _e === void 0 ? void 0 : _e.length) !== null && _f !== void 0 ? _f : 0);
        if (bottomRight.x <= topLeft.x || bottomRight.y <= topLeft.y) {
            return [];
        }
        const bitmap = [];
        for (let y = topLeft.y; y < bottomRight.y; y++) {
            const row = [];
            for (let x = topLeft.x; x < bottomRight.x; x++) {
                const tileData = (_h = (_g = layer.data) === null || _g === void 0 ? void 0 : _g[y]) === null || _h === void 0 ? void 0 : _h[x];
                if (tileData === undefined || tileData === -1) {
                    row.push(false);
                    continue;
                }
                const tile = (_j = layer.tiles) === null || _j === void 0 ? void 0 : _j[tileData];
                if (!tile) {
                    row.push(false);
                    continue;
                }
                if (fieldName && !tile[fieldName]) {
                    row.push(false);
                    continue;
                }
                row.push(true);
            }
            bitmap.push(row);
        }
        return (0, bitmap_decompose_1.bitmapToRectangles)(bitmap);
    }
    /**
     * Get the tile at a given position and in the specified layer
     *
     * If no layer is specified, return a dictionary of layer names to tile
     * definitions (i.e. return all layers)
     *
     * If no tile exists at this position, return null
     */
    getTileAtPosition(position, layerName) {
        if (layerName) {
            return this.getTileAtPositionInLayer(position, layerName);
        }
        const result = {};
        for (const layer of this.options.layers) {
            result[layer.name] = this.getTileAtPositionInLayer(position, layer.name);
        }
        return result;
    }
    getTileAtPositionInLayer(position, layerName) {
        var _a, _b, _c;
        const tilePosition = vec_1.vec.map(vec_1.vec.mul(position, 1 / this.options.tileSize), Math.floor);
        const layer = this.options.layers.find((l) => l.name === layerName);
        if (!layer) {
            return null;
        }
        const tileData = (_b = (_a = layer.data) === null || _a === void 0 ? void 0 : _a[tilePosition.y]) === null || _b === void 0 ? void 0 : _b[tilePosition.x];
        if (tileData === undefined || tileData === -1) {
            return null;
        }
        if (layer.tiles) {
            return (_c = layer.tiles[tileData]) !== null && _c !== void 0 ? _c : null;
        }
        return null;
    }
    hashVector(v) {
        return vec_1.vec.str(v);
    }
    draw(context, screen, position, scale) {
        var _a, _b, _c, _d;
        const absoluteChunkSize = this.options.tileSize * this.options.chunkSize;
        const chunkBorder = (0, vec_1.vec)(this.options.chunkBorder);
        // Maybe clamp scale
        let actualScale = scale;
        if (this.options.minScale && actualScale < this.options.minScale) {
            actualScale = this.options.minScale;
        }
        if (this.options.maxScale && actualScale > this.options.maxScale) {
            actualScale = this.options.maxScale;
        }
        // Maybe clamp position to bounds
        let actualPosition = (0, vec_1.vec)(position);
        if (this.options.bounds && this.options.clampPositionToBounds) {
            const tileSizeScaled = this.options.tileSize / actualScale;
            const halfScreenScaled = vec_1.vec.map(vec_1.vec.mul(screen, 1 / (actualScale * 2)), Math.ceil);
            const minPosition = (0, vec_1.vec)(this.options.bounds.topLeft.x * tileSizeScaled + halfScreenScaled.x, this.options.bounds.topLeft.y * tileSizeScaled + halfScreenScaled.y);
            const maxPosition = (0, vec_1.vec)(this.options.bounds.bottomRight.x * tileSizeScaled - halfScreenScaled.x, this.options.bounds.bottomRight.y * tileSizeScaled - halfScreenScaled.y);
            actualPosition = (0, vec_1.vec)((0, utils_1.clamp)(actualPosition.x, minPosition.x, maxPosition.x), (0, utils_1.clamp)(actualPosition.y, minPosition.y, maxPosition.y));
        }
        const screenSizeInChunks = vec_1.vec.map(vec_1.vec.mul(screen, 1 / (absoluteChunkSize * actualScale)), Math.ceil);
        const screenCenterChunk = vec_1.vec.map(vec_1.vec.mul(actualPosition, 1 / absoluteChunkSize), Math.floor);
        const topLeftChunk = vec_1.vec.sub(vec_1.vec.sub(screenCenterChunk, vec_1.vec.map(vec_1.vec.mul(screenSizeInChunks, 0.5), Math.ceil)), chunkBorder);
        const bottomRightChunk = vec_1.vec.add(vec_1.vec.add(screenCenterChunk, vec_1.vec.map(vec_1.vec.mul(screenSizeInChunks, 0.5), Math.ceil)), chunkBorder);
        context.save();
        context.scale(actualScale, actualScale);
        context.translate(-actualPosition.x + screen.x / (actualScale * 2), -actualPosition.y + screen.y / (actualScale * 2));
        (_b = (_a = this.options).preRender) === null || _b === void 0 ? void 0 : _b.call(_a, context, this, screen, actualPosition, actualScale);
        // Render chunks
        for (let y = topLeftChunk.y; y < bottomRightChunk.y; y++) {
            for (let x = topLeftChunk.x; x < bottomRightChunk.x; x++) {
                const chunkPosition = (0, vec_1.vec)(x, y);
                const chunkAbsolutePosition = vec_1.vec.mul(chunkPosition, absoluteChunkSize);
                // Check if we have this chunk in the cache
                const chunkHash = this.hashVector(chunkPosition);
                if (!this.chunkBuffer.has(chunkHash)) {
                    this.chunkBuffer.set(chunkHash, this.generateChunk(chunkPosition, absoluteChunkSize));
                }
                const chunk = this.chunkBuffer.get(chunkHash);
                if (chunk) {
                    context.drawImage(chunk.image, chunkAbsolutePosition.x, chunkAbsolutePosition.y);
                }
            }
        }
        (_d = (_c = this.options).postRender) === null || _d === void 0 ? void 0 : _d.call(_c, context, this, screen, actualPosition, actualScale);
        // Render debug helpers
        if (this.options.debug.showTileBorders) {
            const topLeftTile = vec_1.vec.mul(vec_1.vec.sub(screenCenterChunk, vec_1.vec.add(vec_1.vec.map(vec_1.vec.mul(screenSizeInChunks, 0.5), Math.ceil), (0, vec_1.vec)(1))), this.options.chunkSize);
            const bottomRightTile = vec_1.vec.mul(vec_1.vec.add(screenCenterChunk, vec_1.vec.add(vec_1.vec.map(vec_1.vec.mul(screenSizeInChunks, 0.5), Math.ceil), (0, vec_1.vec)(1))), this.options.chunkSize);
            for (let y = topLeftTile.y; y < bottomRightTile.y; y++) {
                this.drawLine(context, (0, vec_1.vec)(actualPosition.x - screen.x / (actualScale * 2), y * this.options.tileSize), (0, vec_1.vec)(actualPosition.x + screen.x / (actualScale * 2), y * this.options.tileSize), TileMap.DEBUG_TILE_BORDER_COLOUR, TileMap.DEBUG_TILE_BORDER_LINE_WIDTH);
            }
            for (let x = topLeftTile.x; x < bottomRightTile.x; x++) {
                this.drawLine(context, (0, vec_1.vec)(x * this.options.tileSize, actualPosition.y - screen.y / (actualScale * 2)), (0, vec_1.vec)(x * this.options.tileSize, actualPosition.y + screen.y / (actualScale * 2)), TileMap.DEBUG_TILE_BORDER_COLOUR, TileMap.DEBUG_TILE_BORDER_LINE_WIDTH);
            }
        }
        if (this.options.debug.showChunkBorders) {
            for (let y = topLeftChunk.y; y < bottomRightChunk.y; y++) {
                this.drawLine(context, (0, vec_1.vec)(actualPosition.x - screen.x / (actualScale * 2), y * absoluteChunkSize), (0, vec_1.vec)(actualPosition.x + screen.x / (actualScale * 2), y * absoluteChunkSize), TileMap.DEBUG_CHUNK_BORDER_COLOUR, TileMap.DEBUG_CHUNK_BORDER_LINE_WIDTH);
            }
            for (let x = topLeftChunk.x; x < bottomRightChunk.x; x++) {
                this.drawLine(context, (0, vec_1.vec)(x * absoluteChunkSize, actualPosition.y - screen.y / (actualScale * 2)), (0, vec_1.vec)(x * absoluteChunkSize, actualPosition.y + screen.y / (actualScale * 2)), TileMap.DEBUG_CHUNK_BORDER_COLOUR, TileMap.DEBUG_CHUNK_BORDER_LINE_WIDTH);
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
                    context.fillText(`${x}, ${y}`, x * absoluteChunkSize + absoluteChunkSize / 2, y * absoluteChunkSize + absoluteChunkSize / 2);
                }
            }
            context.restore();
        }
        if (this.options.debug.showOrigin &&
            pointInRectangle((0, vec_1.vec)(0, 0), topLeftChunk, bottomRightChunk)) {
            this.drawCross(context, (0, vec_1.vec)(0, 0), TileMap.DEBUG_ORIGIN_COLOUR, TileMap.DEBUG_ORIGIN_LINE_WIDTH, TileMap.DEBUG_ORIGIN_SIZE);
        }
        context.restore();
    }
    generateChunk(chunkPosition, absoluteChunkSize) {
        var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k, _l, _m;
        const chunkCanvas = document.createElement('canvas');
        const chunkContext = chunkCanvas.getContext('2d');
        chunkCanvas.width = absoluteChunkSize;
        chunkCanvas.height = absoluteChunkSize;
        let chunk = {
            chunkPosition,
            image: chunkCanvas,
        };
        const topLeftTile = vec_1.vec.mul(chunkPosition, this.options.chunkSize);
        const bottomRightTile = vec_1.vec.add(topLeftTile, (0, vec_1.vec)(this.options.chunkSize - 1));
        const boundsTopLeft = (_b = (_a = this.options.bounds) === null || _a === void 0 ? void 0 : _a.topLeft) !== null && _b !== void 0 ? _b : (0, vec_1.vec)(0);
        if (this.options.preGenerateChunk) {
            const result = this.options.preGenerateChunk(chunkContext, this, {
                topLeft: topLeftTile,
                bottomRight: bottomRightTile,
            }, chunkPosition);
            if (Array.isArray(result)) {
                if (!result[1]) {
                    return chunk;
                }
            }
        }
        // Default generation, render tiles from tilemap data
        for (const layer of this.options.layers) {
            chunkContext.save();
            chunkContext.globalAlpha = (_c = layer.opacity) !== null && _c !== void 0 ? _c : 1;
            const alignment = (_d = layer.alignment) !== null && _d !== void 0 ? _d : TileAlignment.Center;
            for (let y = topLeftTile.y; y <= bottomRightTile.y; y++) {
                for (let x = topLeftTile.x; x <= bottomRightTile.x; x++) {
                    const tilePosition = (0, vec_1.vec)(x, y);
                    (_e = layer.preRenderTile) === null || _e === void 0 ? void 0 : _e.call(layer, chunkContext, this, layer, chunkPosition, tilePosition);
                    const tileDataPosition = vec_1.vec.sub(tilePosition, boundsTopLeft);
                    if (tileDataPosition.x < 0 || tileDataPosition.y < 0) {
                        continue;
                    }
                    const tileData = (_g = (_f = layer.data) === null || _f === void 0 ? void 0 : _f[tileDataPosition.y]) === null || _g === void 0 ? void 0 : _g[tileDataPosition.x];
                    if (tileData === undefined || tileData === -1) {
                        continue;
                    }
                    const tileImage = (_j = (_h = layer.tiles) === null || _h === void 0 ? void 0 : _h[tileData]) === null || _j === void 0 ? void 0 : _j.image;
                    if (!tileImage) {
                        continue;
                    }
                    const tileAbsolutePosition = vec_1.vec.sub(vec_1.vec.mul(tilePosition, this.options.tileSize), vec_1.vec.mul(chunkPosition, absoluteChunkSize));
                    // Tile clipping
                    if (layer.clip) {
                        chunkContext.save();
                        chunkContext.beginPath();
                        chunkContext.rect(tileAbsolutePosition.x, tileAbsolutePosition.y, this.options.tileSize, this.options.tileSize);
                        chunkContext.clip();
                    }
                    // Tile alignment
                    let tileImageAbsolutePosition;
                    switch (alignment) {
                        case TileAlignment.TopLeft:
                            tileImageAbsolutePosition = (0, vec_1.vec)(tileAbsolutePosition);
                            break;
                        case TileAlignment.Top:
                            tileImageAbsolutePosition = (0, vec_1.vec)((tileAbsolutePosition.x + this.options.tileSize / 2) - tileImage.width / 2, tileAbsolutePosition.y);
                            break;
                        case TileAlignment.TopRight:
                            tileImageAbsolutePosition = (0, vec_1.vec)(tileAbsolutePosition.x + this.options.tileSize - tileImage.width, tileAbsolutePosition.y);
                            break;
                        case TileAlignment.Left:
                            tileImageAbsolutePosition = (0, vec_1.vec)(tileAbsolutePosition.x, (tileAbsolutePosition.y + this.options.tileSize / 2) - tileImage.height / 2);
                            break;
                        case TileAlignment.Center:
                            tileImageAbsolutePosition = (0, vec_1.vec)((tileAbsolutePosition.x + this.options.tileSize / 2) - tileImage.width / 2, (tileAbsolutePosition.y + this.options.tileSize / 2) - tileImage.height / 2);
                            break;
                        case TileAlignment.Right:
                            tileImageAbsolutePosition = (0, vec_1.vec)(tileAbsolutePosition.x + this.options.tileSize - tileImage.width, (tileAbsolutePosition.y + this.options.tileSize / 2) - tileImage.height / 2);
                            break;
                        case TileAlignment.BottomLeft:
                            tileImageAbsolutePosition = (0, vec_1.vec)(tileAbsolutePosition.x, tileAbsolutePosition.y + this.options.tileSize - tileImage.height);
                            break;
                        case TileAlignment.Bottom:
                            tileImageAbsolutePosition = (0, vec_1.vec)((tileAbsolutePosition.x + this.options.tileSize / 2) - tileImage.width / 2, tileAbsolutePosition.y + this.options.tileSize - tileImage.height);
                            break;
                        case TileAlignment.BottomRight:
                            tileImageAbsolutePosition = (0, vec_1.vec)(tileAbsolutePosition.x + this.options.tileSize - tileImage.width, tileAbsolutePosition.y + this.options.tileSize - tileImage.height);
                            break;
                    }
                    chunkContext.drawImage(tileImage, tileImageAbsolutePosition.x, tileImageAbsolutePosition.y);
                    if (layer.clip) {
                        chunkContext.restore();
                    }
                    (_k = layer.postRenderTile) === null || _k === void 0 ? void 0 : _k.call(layer, chunkCanvas, chunkContext, this, layer, chunkPosition, tilePosition);
                }
            }
            chunkContext.restore();
        }
        (_m = (_l = this.options).postGenerateChunk) === null || _m === void 0 ? void 0 : _m.call(_l, chunkCanvas, chunkContext, this, {
            topLeft: topLeftTile,
            bottomRight: bottomRightTile,
        }, chunkPosition);
        return chunk;
    }
    drawLine(context, start, end, colour, lineWidth) {
        context.save();
        context.lineWidth = lineWidth;
        context.strokeStyle = colour;
        context.beginPath();
        context.moveTo(start.x, start.y);
        context.lineTo(end.x, end.y);
        context.stroke();
        context.restore();
    }
    drawCross(context, position, colour, lineWidth, size) {
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
exports.TileMap = TileMap;
TileMap.DEFAULT_OPTIONS = {
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
TileMap.DEBUG_ORIGIN_COLOUR = 'cyan';
TileMap.DEBUG_ORIGIN_LINE_WIDTH = 2;
TileMap.DEBUG_ORIGIN_SIZE = 10;
TileMap.DEBUG_CHUNK_BORDER_COLOUR = 'yellow';
TileMap.DEBUG_CHUNK_BORDER_LINE_WIDTH = 2;
TileMap.DEBUG_CHUNK_LABEL_COLOUR = 'white';
TileMap.DEBUG_CHUNK_LABEL_FONT = '12px monospace';
TileMap.DEBUG_TILE_BORDER_COLOUR = 'orange';
TileMap.DEBUG_TILE_BORDER_LINE_WIDTH = 1;
/**
 * Content Manager Processor wrapper which converts TileMapOptionsData into
 * TileMapOptions
 *
 * @see https://www.npmjs.com/package/@basementuniverse/content-manager
 */
async function tileMapOptionsContentProcessor(content, data, options) {
    const getImageFromContent = (name) => {
        var _a;
        const image = (_a = content[name]) === null || _a === void 0 ? void 0 : _a.content;
        if (!image) {
            throw new Error(`Image '${name}' not found`);
        }
        return image;
    };
    const result = data;
    if (result.layers) {
        for (const [i, layer] of result.layers.entries()) {
            // Replace imageName with image in the tile definitions array
            if (layer.tiles) {
                for (const [j, tile] of layer.tiles.entries()) {
                    result.layers[i].tiles[j].image = getImageFromContent(tile.imageName);
                    delete result.layers[i].tiles[j].imageName;
                }
            }
            // Decompress layer data
            if ((options === null || options === void 0 ? void 0 : options.decompressData) && layer.data && layer.width) {
                result.layers[i].data = (0, utils_1.chunk)((0, decode_1.default)(layer.data), layer.width);
                delete result.layers[i].width;
            }
        }
    }
    // @ts-ignore
    data.content = result;
}
exports.tileMapOptionsContentProcessor = tileMapOptionsContentProcessor;
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiaW5kZXguanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyIuLi9pbmRleC50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7Ozs7QUFBQSxxQ0FBaUM7QUFDakMsNkRBQXFDO0FBQ3JDLCtDQUE0QztBQUM1QyxtREFBdUQ7QUFDdkQseURBQW1FO0FBNFFuRSxJQUFZLGFBVVg7QUFWRCxXQUFZLGFBQWE7SUFDdkIsdURBQVcsQ0FBQTtJQUNYLCtDQUFHLENBQUE7SUFDSCx5REFBUSxDQUFBO0lBQ1IsaURBQUksQ0FBQTtJQUNKLHFEQUFNLENBQUE7SUFDTixtREFBSyxDQUFBO0lBQ0wsNkRBQVUsQ0FBQTtJQUNWLHFEQUFNLENBQUE7SUFDTiwrREFBVyxDQUFBO0FBQ2IsQ0FBQyxFQVZXLGFBQWEsR0FBYixxQkFBYSxLQUFiLHFCQUFhLFFBVXhCO0FBZUQsU0FBUyxnQkFBZ0IsQ0FDdkIsS0FBVSxFQUNWLE9BQVksRUFDWixXQUFnQjtJQUVoQixPQUFPLENBQ0wsS0FBSyxDQUFDLENBQUMsSUFBSSxPQUFPLENBQUMsQ0FBQztRQUNwQixLQUFLLENBQUMsQ0FBQyxJQUFJLE9BQU8sQ0FBQyxDQUFDO1FBQ3BCLEtBQUssQ0FBQyxDQUFDLEdBQUcsV0FBVyxDQUFDLENBQUM7UUFDdkIsS0FBSyxDQUFDLENBQUMsR0FBRyxXQUFXLENBQUMsQ0FBQyxDQUN4QixDQUFDO0FBQ0osQ0FBQztBQUVELE1BQWEsT0FBTztJQWlDbEIsWUFBbUIsT0FBb0M7UUFDckQsTUFBTSxhQUFhLEdBQUcsTUFBTSxDQUFDLE1BQU0sQ0FDakMsRUFBRSxFQUNGLE9BQU8sQ0FBQyxlQUFlLEVBQ3ZCLE9BQU8sYUFBUCxPQUFPLGNBQVAsT0FBTyxHQUFJLEVBQUUsQ0FDZCxDQUFDO1FBRUYsSUFBSSxDQUFDLGFBQWEsQ0FBQyxLQUFLLElBQUksYUFBYSxDQUFDLEtBQUssS0FBSyxJQUFJLEVBQUU7WUFDeEQsYUFBYSxDQUFDLEtBQUssR0FBRztnQkFDcEIsVUFBVSxFQUFFLENBQUMs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

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	
/******/ 	// startup
/******/ 	// Load entry module and return exports
/******/ 	// This entry module is referenced by other modules so it can't be inlined
/******/ 	var __webpack_exports__ = __webpack_require__("./index.ts");
/******/ 	
/******/ 	return __webpack_exports__;
/******/ })()
;
});
//# sourceMappingURL=data:application/json;charset=utf-8;base64,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QTtBQUNBLGVBQWU7QUFDZjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxTQUFTO0FBQ1Q7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxzQ0FBc0MsS0FBSztBQUMzQztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0Esc0NBQXNDO0FBQ3RDLDJDQUEyQzs7Ozs7O1VDdFkzQztVQUNBOztVQUVBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBOztVQUVBO1VBQ0E7O1VBRUE7VUFDQTtVQUNBOzs7O1VFdEJBO1VBQ0E7VUFDQTtVQUNBIiwic291cmNlcyI6WyJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2UvdGlsZS1tYXAvd2VicGFjay91bml2ZXJzYWxNb2R1bGVEZWZpbml0aW9uIiwid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3RpbGUtbWFwLy4vbm9kZV9tb2R1bGVzL0BiYXNlbWVudHVuaXZlcnNlL3V0aWxzL3V0aWxzLmpzIiwid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3RpbGUtbWFwLy4vbm9kZV9tb2R1bGVzL0BiYXNlbWVudHVuaXZlcnNlL3ZlYy92ZWMuanMiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2UvdGlsZS1tYXAvLi9ub2RlX21vZHVsZXMvZmFzdC1ybGUvZGVjb2RlLmpzIiwid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3RpbGUtbWFwLy4vbm9kZV9tb2R1bGVzL2xydV9tYXAvZGlzdC9scnUuanMiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2UvdGlsZS1tYXAvLi9iaXRtYXAtZGVjb21wb3NlLnRzIiwid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3RpbGUtbWFwLy4vaW5kZXgudHMiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2UvdGlsZS1tYXAvd2VicGFjay9ib290c3RyYXAiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2UvdGlsZS1tYXAvd2VicGFjay9iZWZvcmUtc3RhcnR1cCIsIndlYnBhY2s6Ly9AYmFzZW1lbnR1bml2ZXJzZS90aWxlLW1hcC93ZWJwYWNrL3N0YXJ0dXAiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2UvdGlsZS1tYXAvd2VicGFjay9hZnRlci1zdGFydHVwIl0sInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiB3ZWJwYWNrVW5pdmVyc2FsTW9kdWxlRGVmaW5pdGlvbihyb290LCBmYWN0b3J5KSB7XG5cdGlmKHR5cGVvZiBleHBvcnRzID09PSAnb2JqZWN0JyAmJiB0eXBlb2YgbW9kdWxlID09PSAnb2JqZWN0Jylcblx0XHRtb2R1bGUuZXhwb3J0cyA9IGZhY3RvcnkoKTtcblx0ZWxzZSBpZih0eXBlb2YgZGVmaW5lID09PSAnZnVuY3Rpb24nICYmIGRlZmluZS5hbWQpXG5cdFx0ZGVmaW5lKFtdLCBmYWN0b3J5KTtcblx0ZWxzZSB7XG5cdFx0dmFyIGEgPSBmYWN0b3J5KCk7XG5cdFx0Zm9yKHZhciBpIGluIGEpICh0eXBlb2YgZXhwb3J0cyA9PT0gJ29iamVjdCcgPyBleHBvcnRzIDogcm9vdClbaV0gPSBhW2ldO1xuXHR9XG59KShzZWxmLCAoKSA9PiB7XG5yZXR1cm4gIiwiLyoqXG4gKiBAb3ZlcnZpZXcgQSBsaWJyYXJ5IG9mIHVzZWZ1bCBmdW5jdGlvbnNcbiAqIEBhdXRob3IgR29yZG9uIExhcnJpZ2FuXG4gKi9cblxuLyoqXG4gKiBDaGVjayBpZiB0d28gbnVtYmVycyBhcmUgYXBwcm94aW1hdGVseSBlcXVhbFxuICogQHBhcmFtIHtudW1iZXJ9IGEgTnVtYmVyIGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBiIE51bWJlciBiXG4gKiBAcGFyYW0ge251bWJlcn0gW3A9TnVtYmVyLkVQU0lMT05dIFRoZSBwcmVjaXNpb24gdmFsdWVcbiAqIEByZXR1cm4ge2Jvb2xlYW59IFRydWUgaWYgbnVtYmVycyBhIGFuZCBiIGFyZSBhcHByb3hpbWF0ZWx5IGVxdWFsXG4gKi9cbmNvbnN0IGZsb2F0RXF1YWxzID0gKGEsIGIsIHAgPSBOdW1iZXIuRVBTSUxPTikgPT4gTWF0aC5hYnMoYSAtIGIpIDwgcDtcblxuLyoqXG4gKiBDbGFtcCBhIG51bWJlciBiZXR3ZWVuIG1pbiBhbmQgbWF4XG4gKiBAcGFyYW0ge251bWJlcn0gYSBUaGUgbnVtYmVyIHRvIGNsYW1wXG4gKiBAcGFyYW0ge251bWJlcn0gW21pbj0wXSBUaGUgbWluaW11bSB2YWx1ZVxuICogQHBhcmFtIHtudW1iZXJ9IFttYXg9MV0gVGhlIG1heGltdW0gdmFsdWVcbiAqIEByZXR1cm4ge251bWJlcn0gQSBjbGFtcGVkIG51bWJlclxuICovXG5jb25zdCBjbGFtcCA9IChhLCBtaW4gPSAwLCBtYXggPSAxKSA9PiBhIDwgbWluID8gbWluIDogKGEgPiBtYXggPyBtYXggOiBhKTtcblxuLyoqXG4gKiBHZXQgdGhlIGZyYWN0aW9uYWwgcGFydCBvZiBhIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGEgVGhlIG51bWJlciBmcm9tIHdoaWNoIHRvIGdldCB0aGUgZnJhY3Rpb25hbCBwYXJ0XG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSBmcmFjdGlvbmFsIHBhcnQgb2YgdGhlIG51bWJlclxuICovXG5jb25zdCBmcmFjID0gYSA9PiBhID49IDAgPyBhIC0gTWF0aC5mbG9vcihhKSA6IGEgLSBNYXRoLmNlaWwoYSk7XG5cbi8qKlxuICogUm91bmQgbiB0byBkIGRlY2ltYWwgcGxhY2VzXG4gKiBAcGFyYW0ge251bWJlcn0gbiBUaGUgbnVtYmVyIHRvIHJvdW5kXG4gKiBAcGFyYW0ge251bWJlcn0gW2Q9MF0gVGhlIG51bWJlciBvZiBkZWNpbWFsIHBsYWNlcyB0byByb3VuZCB0b1xuICogQHJldHVybiB7bnVtYmVyfSBBIHJvdW5kZWQgbnVtYmVyXG4gKi9cbmNvbnN0IHJvdW5kID0gKG4sIGQgPSAwKSA9PiB7XG4gIGNvbnN0IHAgPSBNYXRoLnBvdygxMCwgZCk7XG4gIHJldHVybiBNYXRoLnJvdW5kKG4gKiBwICsgTnVtYmVyLkVQU0lMT04pIC8gcDtcbn1cblxuLyoqXG4gKiBEbyBhIGxpbmVhciBpbnRlcnBvbGF0aW9uIGJldHdlZW4gYSBhbmQgYlxuICogQHBhcmFtIHtudW1iZXJ9IGEgVGhlIG1pbmltdW0gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gYiBUaGUgbWF4aW11bSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSBpbnRlcnBvbGF0aW9uIHZhbHVlLCBzaG91bGQgYmUgaW4gdGhlIGludGVydmFsIFswLCAxXVxuICogQHJldHVybiB7bnVtYmVyfSBBbiBpbnRlcnBvbGF0ZWQgdmFsdWUgaW4gdGhlIGludGVydmFsIFthLCBiXVxuICovXG5jb25zdCBsZXJwID0gKGEsIGIsIGkpID0+IGEgKyAoYiAtIGEpICogaTtcblxuLyoqXG4gKiBHZXQgdGhlIHBvc2l0aW9uIG9mIGkgYmV0d2VlbiBhIGFuZCBiXG4gKiBAcGFyYW0ge251bWJlcn0gYSBUaGUgbWluaW11bSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBiIFRoZSBtYXhpbXVtIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIGludGVycG9sYXRlZCB2YWx1ZSBpbiB0aGUgaW50ZXJ2YWwgW2EsIGJdXG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSBwb3NpdGlvbiBvZiBpIGJldHdlZW4gYSBhbmQgYlxuICovXG5jb25zdCB1bmxlcnAgPSAoYSwgYiwgaSkgPT4gKGkgLSBhKSAvIChiIC0gYSk7XG5cbi8qKlxuICogRG8gYSBiaWxpbmVhciBpbnRlcnBvbGF0aW9uXG4gKiBAcGFyYW0ge251bWJlcn0gYzAwIFRvcC1sZWZ0IHZhbHVlXG4gKiBAcGFyYW0ge251bWJlcn0gYzEwIFRvcC1yaWdodCB2YWx1ZVxuICogQHBhcmFtIHtudW1iZXJ9IGMwMSBCb3R0b20tbGVmdCB2YWx1ZVxuICogQHBhcmFtIHtudW1iZXJ9IGMxMSBCb3R0b20tcmlnaHQgdmFsdWVcbiAqIEBwYXJhbSB7bnVtYmVyfSBpeCBJbnRlcnBvbGF0aW9uIHZhbHVlIGFsb25nIHhcbiAqIEBwYXJhbSB7bnVtYmVyfSBpeSBJbnRlcnBvbGF0aW9uIHZhbHVlIGFsb25nIHlcbiAqIEByZXR1cm4ge251bWJlcn0gQSBiaWxpbmVhciBpbnRlcnBvbGF0ZWQgdmFsdWVcbiAqL1xuY29uc3QgYmxlcnAgPSAoYzAwLCBjMTAsIGMwMSwgYzExLCBpeCwgaXkpID0+IGxlcnAobGVycChjMDAsIGMxMCwgaXgpLCBsZXJwKGMwMSwgYzExLCBpeCksIGl5KTtcblxuLyoqXG4gKiBSZS1tYXAgYSBudW1iZXIgaSBmcm9tIHJhbmdlIGExLi4uYTIgdG8gYjEuLi5iMlxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIG51bWJlciB0byByZS1tYXBcbiAqIEBwYXJhbSB7bnVtYmVyfSBhMVxuICogQHBhcmFtIHtudW1iZXJ9IGEyXG4gKiBAcGFyYW0ge251bWJlcn0gYjFcbiAqIEBwYXJhbSB7bnVtYmVyfSBiMlxuICogQHJldHVybiB7bnVtYmVyfVxuICovXG5jb25zdCByZW1hcCA9IChpLCBhMSwgYTIsIGIxLCBiMikgPT4gYjEgKyAoaSAtIGExKSAqIChiMiAtIGIxKSAvIChhMiAtIGExKTtcblxuLyoqXG4gKiBEbyBhIHNtb290aCBpbnRlcnBvbGF0aW9uIGJldHdlZW4gYSBhbmQgYlxuICogQHBhcmFtIHtudW1iZXJ9IGEgVGhlIG1pbmltdW0gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gYiBUaGUgbWF4aW11bSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSBpbnRlcnBvbGF0aW9uIHZhbHVlXG4gKiBAcmV0dXJuIHtudW1iZXJ9IEFuIGludGVycG9sYXRlZCB2YWx1ZSBpbiB0aGUgaW50ZXJ2YWwgW2EsIGJdXG4gKi9cbmNvbnN0IHNtb290aHN0ZXAgPSAoYSwgYiwgaSkgPT4gbGVycChhLCBiLCAzICogTWF0aC5wb3coaSwgMikgLSAyICogTWF0aC5wb3coaSwgMykpO1xuXG4vKipcbiAqIEdldCBhbiBhbmdsZSBpbiByYWRpYW5zXG4gKiBAcGFyYW0ge251bWJlcn0gZGVncmVlcyBUaGUgYW5nbGUgaW4gZGVncmVlc1xuICogQHJldHVybiB7bnVtYmVyfSBUaGUgYW5nbGUgaW4gcmFkaWFuc1xuICovXG5jb25zdCByYWRpYW5zID0gZGVncmVlcyA9PiAoTWF0aC5QSSAvIDE4MCkgKiBkZWdyZWVzO1xuXG4vKipcbiAqIEdldCBhbiBhbmdsZSBpbiBkZWdyZWVzXG4gKiBAcGFyYW0ge251bWJlcn0gcmFkaWFucyBUaGUgYW5nbGUgaW4gcmFkaWFuc1xuICogQHJldHVybiB7bnVtYmVyfSBUaGUgYW5nbGUgaW4gZGVncmVlc1xuICovXG5jb25zdCBkZWdyZWVzID0gcmFkaWFucyA9PiAoMTgwIC8gTWF0aC5QSSkgKiByYWRpYW5zO1xuXG4vKipcbiAqIEdldCBhIHJhbmRvbSBmbG9hdCBpbiB0aGUgaW50ZXJ2YWwgW21pbiwgbWF4KVxuICogQHBhcmFtIHtudW1iZXJ9IG1pbiBJbmNsdXNpdmUgbWluXG4gKiBAcGFyYW0ge251bWJlcn0gbWF4IEV4Y2x1c2l2ZSBtYXhcbiAqIEByZXR1cm4ge251bWJlcn0gQSByYW5kb20gZmxvYXQgaW4gdGhlIGludGVydmFsIFttaW4sIG1heClcbiAqL1xuY29uc3QgcmFuZG9tQmV0d2VlbiA9IChtaW4sIG1heCkgPT4gTWF0aC5yYW5kb20oKSAqIChtYXggLSBtaW4pICsgbWluO1xuXG4vKipcbiAqIEdldCBhIHJhbmRvbSBpbnRlZ2VyIGluIHRoZSBpbnRlcnZhbCBbbWluLCBtYXhdXG4gKiBAcGFyYW0ge251bWJlcn0gbWluIEluY2x1c2l2ZSBtaW5cbiAqIEBwYXJhbSB7bnVtYmVyfSBtYXggSW5jbHVzaXZlIG1heFxuICogQHJldHVybiB7bnVtYmVyfSBBIHJhbmRvbSBpbnRlZ2VyIGluIHRoZSBpbnRlcnZhbCBbbWluLCBtYXhdXG4gKi9cbmNvbnN0IHJhbmRvbUludEJldHdlZW4gPSAobWluLCBtYXgpID0+IE1hdGguZmxvb3IoTWF0aC5yYW5kb20oKSAqIChtYXggLSBtaW4gKyAxKSkgKyBtaW47XG5cbi8qKlxuICogR2V0IGEgbm9ybWFsbHktZGlzdHJpYnV0ZWQgcmFuZG9tIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IFttdT0wLjVdIFRoZSBtZWFuIHZhbHVlXG4gKiBAcGFyYW0ge251bWJlcn0gW3NpZ21hPTAuNV0gVGhlIHN0YW5kYXJkIGRldmlhdGlvblxuICogQHBhcmFtIHtudW1iZXJ9IFtzYW1wbGVzPTJdIFRoZSBudW1iZXIgb2Ygc2FtcGxlc1xuICogQHJldHVybiB7bnVtYmVyfSBBIG5vcm1hbGx5LWRpc3RyaWJ1dGVkIHJhbmRvbSBudW1iZXJcbiAqL1xuY29uc3QgY2x0UmFuZG9tID0gKG11ID0gMC41LCBzaWdtYSA9IDAuNSwgc2FtcGxlcyA9IDIpID0+IHtcbiAgbGV0IHRvdGFsID0gMDtcbiAgZm9yIChsZXQgaSA9IHNhbXBsZXM7IGktLTspIHtcbiAgICB0b3RhbCArPSBNYXRoLnJhbmRvbSgpO1xuICB9XG4gIHJldHVybiBtdSArICh0b3RhbCAtIHNhbXBsZXMgLyAyKSAvIChzYW1wbGVzIC8gMikgKiBzaWdtYTtcbn07XG5cbi8qKlxuICogR2V0IGEgbm9ybWFsbHktZGlzdHJpYnV0ZWQgcmFuZG9tIGludGVnZXIgaW4gdGhlIGludGVydmFsIFttaW4sIG1heF1cbiAqIEBwYXJhbSB7bnVtYmVyfSBtaW4gSW5jbHVzaXZlIG1pblxuICogQHBhcmFtIHtudW1iZXJ9IG1heCBJbmNsdXNpdmUgbWF4XG4gKiBAcmV0dXJuIHtudW1iZXJ9IEEgbm9ybWFsbHktZGlzdHJpYnV0ZWQgcmFuZG9tIGludGVnZXJcbiAqL1xuY29uc3QgY2x0UmFuZG9tSW50ID0gKG1pbiwgbWF4KSA9PiBNYXRoLmZsb29yKG1pbiArIGNsdFJhbmRvbSgwLjUsIDAuNSwgMikgKiAobWF4ICsgMSAtIG1pbikpO1xuXG4vKipcbiAqIFJldHVybiBhIHdlaWdodGVkIHJhbmRvbSBpbnRlZ2VyXG4gKiBAcGFyYW0ge0FycmF5PG51bWJlcj59IHcgQW4gYXJyYXkgb2Ygd2VpZ2h0c1xuICogQHJldHVybiB7bnVtYmVyfSBBbiBpbmRleCBmcm9tIHdcbiAqL1xuY29uc3Qgd2VpZ2h0ZWRSYW5kb20gPSB3ID0+IHtcbiAgbGV0IHRvdGFsID0gdy5yZWR1Y2UoKGEsIGkpID0+IGEgKyBpLCAwKSwgbiA9IDA7XG4gIGNvbnN0IHIgPSBNYXRoLnJhbmRvbSgpICogdG90YWw7XG4gIHdoaWxlICh0b3RhbCA+IHIpIHtcbiAgICB0b3RhbCAtPSB3W24rK107XG4gIH1cbiAgcmV0dXJuIG4gLSAxO1xufTtcblxuLyoqXG4gKiBBbiBpbnRlcnBvbGF0aW9uIGZ1bmN0aW9uXG4gKiBAY2FsbGJhY2sgSW50ZXJwb2xhdGlvbkZ1bmN0aW9uXG4gKiBAcGFyYW0ge251bWJlcn0gYSBUaGUgbWluaW11bSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBiIFRoZSBtYXhpbXVtIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIGludGVycG9sYXRpb24gdmFsdWUsIHNob3VsZCBiZSBpbiB0aGUgaW50ZXJ2YWwgWzAsIDFdXG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSBpbnRlcnBvbGF0ZWQgdmFsdWUgaW4gdGhlIGludGVydmFsIFthLCBiXVxuICovXG5cbi8qKlxuICogUmV0dXJuIGFuIGludGVycG9sYXRlZCB2YWx1ZSBmcm9tIGFuIGFycmF5XG4gKiBAcGFyYW0ge0FycmF5PG51bWJlcj59IGEgQW4gYXJyYXkgb2YgdmFsdWVzIGludGVycG9sYXRlXG4gKiBAcGFyYW0ge251bWJlcn0gaSBBIG51bWJlciBpbiB0aGUgaW50ZXJ2YWwgWzAsIDFdXG4gKiBAcGFyYW0ge0ludGVycG9sYXRpb25GdW5jdGlvbn0gW2Y9TWF0aC5sZXJwXSBUaGUgaW50ZXJwb2xhdGlvbiBmdW5jdGlvbiB0byB1c2VcbiAqIEByZXR1cm4ge251bWJlcn0gQW4gaW50ZXJwb2xhdGVkIHZhbHVlIGluIHRoZSBpbnRlcnZhbCBbbWluKGEpLCBtYXgoYSldXG4gKi9cbmNvbnN0IGxlcnBBcnJheSA9IChhLCBpLCBmID0gbGVycCkgPT4ge1xuICBjb25zdCBzID0gaSAqIChhLmxlbmd0aCAtIDEpO1xuICBjb25zdCBwID0gY2xhbXAoTWF0aC50cnVuYyhzKSwgMCwgYS5sZW5ndGggLSAxKTtcbiAgcmV0dXJuIGYoYVtwXSB8fCAwLCBhW3AgKyAxXSB8fCAwLCBmcmFjKHMpKTtcbn07XG5cbi8qKlxuICogR2V0IHRoZSBkb3QgcHJvZHVjdCBvZiB0d28gdmVjdG9yc1xuICogQHBhcmFtIHtBcnJheTxudW1iZXI+fSBhIFZlY3RvciBhXG4gKiBAcGFyYW0ge0FycmF5PG51bWJlcj59IGIgVmVjdG9yIGJcbiAqIEByZXR1cm4ge251bWJlcn0gYSDiiJkgYlxuICovXG5jb25zdCBkb3QgPSAoYSwgYikgPT4gYS5yZWR1Y2UoKG4sIHYsIGkpID0+IG4gKyB2ICogYltpXSwgMCk7XG5cbi8qKlxuICogR2V0IHRoZSBmYWN0b3JpYWwgb2YgYSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBhXG4gKiBAcmV0dXJuIHtudW1iZXJ9IGEhXG4gKi9cbmNvbnN0IGZhY3RvcmlhbCA9IGEgPT4ge1xuICBsZXQgcmVzdWx0ID0gMTtcbiAgZm9yIChsZXQgaSA9IDI7IGkgPD0gYTsgaSsrKSB7XG4gICAgcmVzdWx0ICo9IGk7XG4gIH1cbiAgcmV0dXJuIHJlc3VsdDtcbn07XG5cbi8qKlxuICogR2V0IHRoZSBudW1iZXIgb2YgcGVybXV0YXRpb25zIG9mIHIgZWxlbWVudHMgZnJvbSBhIHNldCBvZiBuIGVsZW1lbnRzXG4gKiBAcGFyYW0ge251bWJlcn0gblxuICogQHBhcmFtIHtudW1iZXJ9IHJcbiAqIEByZXR1cm4ge251bWJlcn0gblByXG4gKi9cbmNvbnN0IG5wciA9IChuLCByKSA9PiBmYWN0b3JpYWwobikgLyBmYWN0b3JpYWwobiAtIHIpO1xuXG4vKipcbiAqIEdldCB0aGUgbnVtYmVyIG9mIGNvbWJpbmF0aW9ucyBvZiByIGVsZW1lbnRzIGZyb20gYSBzZXQgb2YgbiBlbGVtZW50c1xuICogQHBhcmFtIHtudW1iZXJ9IG5cbiAqIEBwYXJhbSB7bnVtYmVyfSByXG4gKiBAcmV0dXJuIHtudW1iZXJ9IG5DclxuICovXG5jb25zdCBuY3IgPSAobiwgcikgPT4gZmFjdG9yaWFsKG4pIC8gKGZhY3RvcmlhbChyKSAqIGZhY3RvcmlhbChuIC0gcikpO1xuXG4vKipcbiAqIEdlbmVyYXRlIGFsbCBjb21iaW5hdGlvbnMgb2YgciBlbGVtZW50cyBmcm9tIGFuIGFycmF5XG4gKlxuICogQGV4YW1wbGVcbiAqIGBgYGpzXG4gKiBjb21iaW5hdGlvbnMoWzEsIDIsIDNdLCAyKTtcbiAqIGBgYFxuICpcbiAqIE91dHB1dDpcbiAqIGBgYGpzb25cbiAqIFtcbiAqICAgWzEsIDJdLFxuICogICBbMSwgM10sXG4gKiAgIFsyLCAzXVxuICogXVxuICogYGBgXG4gKiBAcGFyYW0ge0FycmF5PCo+fSBhXG4gKiBAcGFyYW0ge251bWJlcn0gciBUaGUgbnVtYmVyIG9mIGVsZW1lbnRzIHRvIGNob29zZSBpbiBlYWNoIGNvbWJpbmF0aW9uXG4gKiBAcmV0dXJuIHtBcnJheTxBcnJheTwqPj59IEFuIGFycmF5IG9mIGNvbWJpbmF0aW9uIGFycmF5c1xuICovXG5jb25zdCBjb21iaW5hdGlvbnMgPSAoYSwgcikgPT4ge1xuICBpZiAociA9PT0gMSkge1xuICAgIHJldHVybiBhLm1hcChpdGVtID0+IFtpdGVtXSk7XG4gIH1cblxuICByZXR1cm4gYS5yZWR1Y2UoXG4gICAgKGFjYywgaXRlbSwgaSkgPT4gW1xuICAgICAgLi4uYWNjLFxuICAgICAgLi4uY29tYmluYXRpb25zKGEuc2xpY2UoaSArIDEpLCByIC0gMSkubWFwKGMgPT4gW2l0ZW0sIC4uLmNdKSxcbiAgICBdLFxuICAgIFtdXG4gICk7XG59O1xuXG4vKipcbiAqIEdldCBhIGNhcnRlc2lhbiBwcm9kdWN0IG9mIGFycmF5c1xuICpcbiAqIEBleGFtcGxlXG4gKiBgYGBqc1xuICogY2FydGVzaWFuKFsxLCAyLCAzXSwgWydhJywgJ2InXSk7XG4gKiBgYGBcbiAqXG4gKiBPdXRwdXQ6XG4gKiBgYGBqc29uXG4gKiBbXG4gKiAgIFsxLCBcImFcIl0sXG4gKiAgIFsxLCBcImJcIl0sXG4gKiAgIFsyLCBcImFcIl0sXG4gKiAgIFsyLCBcImJcIl0sXG4gKiAgIFszLCBcImFcIl0sXG4gKiAgIFszLCBcImJcIl1cbiAqIF1cbiAqIGBgYFxuICovXG5jb25zdCBjYXJ0ZXNpYW4gPSAoLi4uYXJyKSA9PlxuICBhcnIucmVkdWNlKFxuICAgIChhLCBiKSA9PiBhLmZsYXRNYXAoYyA9PiBiLm1hcChkID0+IFsuLi5jLCBkXSkpLFxuICAgIFtbXV1cbiAgKTtcblxuLyoqXG4gKiBBIGZ1bmN0aW9uIGZvciBnZW5lcmF0aW5nIGFycmF5IHZhbHVlc1xuICogQGNhbGxiYWNrIFRpbWVzRnVuY3Rpb25cbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSBhcnJheSBpbmRleFxuICogQHJldHVybiB7Kn0gVGhlIGFycmF5IHZhbHVlXG4gKi9cblxuLyoqXG4gKiBSZXR1cm4gYSBuZXcgYXJyYXkgd2l0aCBsZW5ndGggbiBieSBjYWxsaW5nIGZ1bmN0aW9uIGYoaSkgb24gZWFjaCBlbGVtZW50XG4gKiBAcGFyYW0ge1RpbWVzRnVuY3Rpb259IGZcbiAqIEBwYXJhbSB7bnVtYmVyfSBuIFRoZSBzaXplIG9mIHRoZSBhcnJheVxuICogQHJldHVybiB7QXJyYXk8Kj59XG4gKi9cbmNvbnN0IHRpbWVzID0gKGYsIG4pID0+IEFycmF5KG4pLmZpbGwoMCkubWFwKChfLCBpKSA9PiBmKGkpKTtcblxuLyoqXG4gKiBSZXR1cm4gYW4gYXJyYXkgY29udGFpbmluZyBudW1iZXJzIDAtPihuIC0gMSlcbiAqIEBwYXJhbSB7bnVtYmVyfSBuIFRoZSBzaXplIG9mIHRoZSBhcnJheVxuICogQHJldHVybiB7QXJyYXk8bnVtYmVyPn0gQW4gYXJyYXkgb2YgaW50ZWdlcnMgMC0+KG4gLSAxKVxuICovXG5jb25zdCByYW5nZSA9IG4gPT4gdGltZXMoaSA9PiBpLCBuKTtcblxuLyoqXG4gKiBaaXAgMiBhcnJheXMgdG9nZXRoZXIsIGkuZS4gKFsxLCAyLCAzXSwgW2EsIGIsIGNdKSA9PiBbWzEsIGFdLCBbMiwgYl0sIFszLCBjXV1cbiAqIEBwYXJhbSB7QXJyYXk8Kj59IGFcbiAqIEBwYXJhbSB7QXJyYXk8Kj59IGJcbiAqIEByZXR1cm4ge0FycmF5PEFycmF5PCo+Pn1cbiAqL1xuY29uc3QgemlwID0gKGEsIGIpID0+IGEubWFwKChrLCBpKSA9PiBbaywgYltpXV0pO1xuXG4vKipcbiAqIFJldHVybiBhcnJheVtpXSB3aXRoIHBvc2l0aXZlIGFuZCBuZWdhdGl2ZSB3cmFwcGluZ1xuICogQHBhcmFtIHtBcnJheTwqPn0gYVxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIHBvc2l0aXZlbHkvbmVnYXRpdmVseSB3cmFwcGVkIGFycmF5IGluZGV4XG4gKiBAcmV0dXJuIHsqfSBBbiBlbGVtZW50IGZyb20gdGhlIGFycmF5XG4gKi9cbmNvbnN0IGF0ID0gKGEsIGkpID0+IGFbaSA8IDAgPyBhLmxlbmd0aCAtIChNYXRoLmFicyhpICsgMSkgJSBhLmxlbmd0aCkgLSAxIDogaSAlIGEubGVuZ3RoXTtcblxuLyoqXG4gKiBSZXR1cm4gdGhlIGxhc3QgZWxlbWVudCBvZiBhbiBhcnJheSB3aXRob3V0IHJlbW92aW5nIGl0XG4gKiBAcGFyYW0ge0FycmF5PCo+fSBhXG4gKiBAcmV0dXJuIHsqfSBUaGUgbGFzdCBlbGVtZW50IGZyb20gdGhlIGFycmF5XG4gKi9cbmNvbnN0IHBlZWsgPSAoYSkgPT4ge1xuICBpZiAoIWEubGVuZ3RoKSB7XG4gICAgcmV0dXJuIHVuZGVmaW5lZDtcbiAgfVxuXG4gIHJldHVybiBhW2EubGVuZ3RoIC0gMV07XG59O1xuXG4vKipcbiAqIENob3AgYW4gYXJyYXkgaW50byBjaHVua3Mgb2Ygc2l6ZSBuXG4gKiBAcGFyYW0ge0FycmF5PCo+fSBhXG4gKiBAcGFyYW0ge251bWJlcn0gbiBUaGUgY2h1bmsgc2l6ZVxuICogQHJldHVybiB7QXJyYXk8QXJyYXk8Kj4+fSBBbiBhcnJheSBvZiBhcnJheSBjaHVua3NcbiAqL1xuY29uc3QgY2h1bmsgPSAoYSwgbikgPT4gdGltZXMoaSA9PiBhLnNsaWNlKGkgKiBuLCBpICogbiArIG4pLCBNYXRoLmNlaWwoYS5sZW5ndGggLyBuKSk7XG5cbi8qKlxuICogUmFuZG9tbHkgc2h1ZmZsZSBhIHNoYWxsb3cgY29weSBvZiBhbiBhcnJheVxuICogQHBhcmFtIHtBcnJheTwqPn0gYVxuICogQHJldHVybiB7QXJyYXk8Kj59IFRoZSBzaHVmZmxlZCBhcnJheVxuICovXG5jb25zdCBzaHVmZmxlID0gYSA9PiBhLnNsaWNlKCkuc29ydCgoKSA9PiBNYXRoLnJhbmRvbSgpIC0gMC41KTtcblxuLyoqXG4gKiBGbGF0dGVuIGFuIG9iamVjdFxuICogQHBhcmFtIHtvYmplY3R9IG9cbiAqIEBwYXJhbSB7c3RyaW5nfSBjb25jYXRlbmF0b3IgVGhlIHN0cmluZyB0byB1c2UgZm9yIGNvbmNhdGVuYXRpbmcga2V5c1xuICogQHJldHVybiB7b2JqZWN0fSBBIGZsYXR0ZW5lZCBvYmplY3RcbiAqL1xuY29uc3QgZmxhdCA9IChvLCBjb25jYXRlbmF0b3IgPSAnLicpID0+IHtcbiAgcmV0dXJuIE9iamVjdC5rZXlzKG8pLnJlZHVjZSgoYWNjLCBrZXkpID0+IHtcbiAgICBpZiAob1trZXldIGluc3RhbmNlb2YgRGF0ZSkge1xuICAgICAgcmV0dXJuIHtcbiAgICAgICAgLi4uYWNjLFxuICAgICAgICBba2V5XTogb1trZXldLnRvSVNPU3RyaW5nKCksXG4gICAgICB9O1xuICAgIH1cblxuICAgIGlmICh0eXBlb2Ygb1trZXldICE9PSAnb2JqZWN0JyB8fCAhb1trZXldKSB7XG4gICAgICByZXR1cm4ge1xuICAgICAgICAuLi5hY2MsXG4gICAgICAgIFtrZXldOiBvW2tleV0sXG4gICAgICB9O1xuICAgIH1cbiAgICBjb25zdCBmbGF0dGVuZWQgPSBmbGF0KG9ba2V5XSwgY29uY2F0ZW5hdG9yKTtcblxuICAgIHJldHVybiB7XG4gICAgICAuLi5hY2MsXG4gICAgICAuLi5PYmplY3Qua2V5cyhmbGF0dGVuZWQpLnJlZHVjZShcbiAgICAgICAgKGNoaWxkQWNjLCBjaGlsZEtleSkgPT4gKHtcbiAgICAgICAgICAuLi5jaGlsZEFjYyxcbiAgICAgICAgICBbYCR7a2V5fSR7Y29uY2F0ZW5hdG9yfSR7Y2hpbGRLZXl9YF06IGZsYXR0ZW5lZFtjaGlsZEtleV0sXG4gICAgICAgIH0pLFxuICAgICAgICB7fVxuICAgICAgKSxcbiAgICB9O1xuICB9LCB7fSk7XG59O1xuXG4vKipcbiAqIFVuZmxhdHRlbiBhbiBvYmplY3RcbiAqIEBwYXJhbSB7b2JqZWN0fSBvXG4gKiBAcGFyYW0ge3N0cmluZ30gY29uY2F0ZW5hdG9yIFRoZSBzdHJpbmcgdG8gY2hlY2sgZm9yIGluIGNvbmNhdGVuYXRlZCBrZXlzXG4gKiBAcmV0dXJuIHtvYmplY3R9IEFuIHVuLWZsYXR0ZW5lZCBvYmplY3RcbiAqL1xuY29uc3QgdW5mbGF0ID0gKG8sIGNvbmNhdGVuYXRvciA9ICcuJykgPT4ge1xuICBsZXQgcmVzdWx0ID0ge30sIHRlbXAsIHN1YnN0cmluZ3MsIHByb3BlcnR5LCBpO1xuXG4gIGZvciAocHJvcGVydHkgaW4gbykge1xuICAgIHN1YnN0cmluZ3MgPSBwcm9wZXJ0eS5zcGxpdChjb25jYXRlbmF0b3IpO1xuICAgIHRlbXAgPSByZXN1bHQ7XG4gICAgZm9yIChpID0gMDsgaSA8IHN1YnN0cmluZ3MubGVuZ3RoIC0gMTsgaSsrKSB7XG4gICAgICBpZiAoIShzdWJzdHJpbmdzW2ldIGluIHRlbXApKSB7XG4gICAgICAgIGlmIChpc0Zpbml0ZShzdWJzdHJpbmdzW2kgKyAxXSkpIHtcbiAgICAgICAgICB0ZW1wW3N1YnN0cmluZ3NbaV1dID0gW107XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgdGVtcFtzdWJzdHJpbmdzW2ldXSA9IHt9O1xuICAgICAgICB9XG4gICAgICB9XG4gICAgICB0ZW1wID0gdGVtcFtzdWJzdHJpbmdzW2ldXTtcbiAgICB9XG4gICAgdGVtcFtzdWJzdHJpbmdzW3N1YnN0cmluZ3MubGVuZ3RoIC0gMV1dID0gb1twcm9wZXJ0eV07XG4gIH1cblxuICByZXR1cm4gcmVzdWx0O1xufTtcblxuLyoqXG4gKiBBIHNwbGl0IHByZWRpY2F0ZVxuICogQGNhbGxiYWNrIFNwbGl0UHJlZGljYXRlXG4gKiBAcGFyYW0ge2FueX0gdmFsdWUgVGhlIGN1cnJlbnQgdmFsdWVcbiAqIEByZXR1cm4ge2Jvb2xlYW59IFRydWUgaWYgdGhlIGFycmF5IHNob3VsZCBzcGxpdCBhdCB0aGlzIGluZGV4XG4gKi9cblxuLyoqXG4gKiBTcGxpdCBhbiBhcnJheSBpbnRvIHN1Yi1hcnJheXMgYmFzZWQgb24gYSBwcmVkaWNhdGVcbiAqIEBwYXJhbSB7QXJyYXk8Kj59IGFycmF5XG4gKiBAcGFyYW0ge1NwbGl0UHJlZGljYXRlfSBwcmVkaWNhdGVcbiAqIEByZXR1cm4ge0FycmF5PEFycmF5PCo+Pn0gQW4gYXJyYXkgb2YgYXJyYXlzXG4gKi9cbmNvbnN0IHNwbGl0ID0gKGFycmF5LCBwcmVkaWNhdGUpID0+IHtcbiAgY29uc3QgcmVzdWx0ID0gW107XG4gIGxldCBjdXJyZW50ID0gW107XG4gIGZvciAoY29uc3QgdmFsdWUgb2YgYXJyYXkpIHtcbiAgICBpZiAocHJlZGljYXRlKHZhbHVlKSkge1xuICAgICAgaWYgKGN1cnJlbnQubGVuZ3RoKSB7XG4gICAgICAgIHJlc3VsdC5wdXNoKGN1cnJlbnQpO1xuICAgICAgfVxuICAgICAgY3VycmVudCA9IFt2YWx1ZV07XG4gICAgfSBlbHNlIHtcbiAgICAgIGN1cnJlbnQucHVzaCh2YWx1ZSk7XG4gICAgfVxuICB9XG4gIHJlc3VsdC5wdXNoKGN1cnJlbnQpO1xuXG4gIHJldHVybiByZXN1bHQ7XG59O1xuXG4vKipcbiAqIFBsdWNrIGtleXMgZnJvbSBhbiBvYmplY3RcbiAqIEBwYXJhbSB7b2JqZWN0fSBvXG4gKiBAcGFyYW0gey4uLnN0cmluZ30ga2V5cyBUaGUga2V5cyB0byBwbHVjayBmcm9tIHRoZSBvYmplY3RcbiAqIEByZXR1cm4ge29iamVjdH0gQW4gb2JqZWN0IGNvbnRhaW5pbmcgdGhlIHBsdWNrZWQga2V5c1xuICovXG5jb25zdCBwbHVjayA9IChvLCAuLi5rZXlzKSA9PiB7XG4gIHJldHVybiBrZXlzLnJlZHVjZShcbiAgICAocmVzdWx0LCBrZXkpID0+IE9iamVjdC5hc3NpZ24ocmVzdWx0LCB7IFtrZXldOiBvW2tleV0gfSksXG4gICAge31cbiAgKTtcbn07XG5cbi8qKlxuICogRXhjbHVkZSBrZXlzIGZyb20gYW4gb2JqZWN0XG4gKiBAcGFyYW0ge29iamVjdH0gb1xuICogQHBhcmFtIHsuLi5zdHJpbmd9IGtleXMgVGhlIGtleXMgdG8gZXhjbHVkZSBmcm9tIHRoZSBvYmplY3RcbiAqIEByZXR1cm4ge29iamVjdH0gQW4gb2JqZWN0IGNvbnRhaW5pbmcgYWxsIGtleXMgZXhjZXB0IGV4Y2x1ZGVkIGtleXNcbiAqL1xuY29uc3QgZXhjbHVkZSA9IChvLCAuLi5rZXlzKSA9PiB7XG4gIHJldHVybiBPYmplY3QuZnJvbUVudHJpZXMoXG4gICAgT2JqZWN0LmVudHJpZXMobykuZmlsdGVyKChba2V5XSkgPT4gIWtleXMuaW5jbHVkZXMoa2V5KSlcbiAgKTtcbn07XG5cbmlmICh0eXBlb2YgbW9kdWxlICE9PSAndW5kZWZpbmVkJykge1xuICBtb2R1bGUuZXhwb3J0cyA9IHtcbiAgICBmbG9hdEVxdWFscyxcbiAgICBjbGFtcCxcbiAgICBmcmFjLFxuICAgIHJvdW5kLFxuICAgIGxlcnAsXG4gICAgdW5sZXJwLFxuICAgIGJsZXJwLFxuICAgIHJlbWFwLFxuICAgIHNtb290aHN0ZXAsXG4gICAgcmFkaWFucyxcbiAgICBkZWdyZWVzLFxuICAgIHJhbmRvbUJldHdlZW4sXG4gICAgcmFuZG9tSW50QmV0d2VlbixcbiAgICBjbHRSYW5kb20sXG4gICAgY2x0UmFuZG9tSW50LFxuICAgIHdlaWdodGVkUmFuZG9tLFxuICAgIGxlcnBBcnJheSxcbiAgICBkb3QsXG4gICAgZmFjdG9yaWFsLFxuICAgIG5wcixcbiAgICBuY3IsXG4gICAgY29tYmluYXRpb25zLFxuICAgIGNhcnRlc2lhbixcbiAgICB0aW1lcyxcbiAgICByYW5nZSxcbiAgICB6aXAsXG4gICAgYXQsXG4gICAgcGVlayxcbiAgICBjaHVuayxcbiAgICBzaHVmZmxlLFxuICAgIGZsYXQsXG4gICAgdW5mbGF0LFxuICAgIHNwbGl0LFxuICAgIHBsdWNrLFxuICAgIGV4Y2x1ZGUsXG4gIH07XG59XG4iLCJjb25zdCB7IHRpbWVzLCBjaHVuaywgZG90IH0gPSByZXF1aXJlKCdAYmFzZW1lbnR1bml2ZXJzZS91dGlscycpO1xuXG4vKipcbiAqIEBvdmVydmlldyBBIHNtYWxsIHZlY3RvciBhbmQgbWF0cml4IGxpYnJhcnlcbiAqIEBhdXRob3IgR29yZG9uIExhcnJpZ2FuXG4gKi9cblxuLyoqXG4gKiBBIDJkIHZlY3RvclxuICogQHR5cGVkZWYge09iamVjdH0gdmVjXG4gKiBAcHJvcGVydHkge251bWJlcn0geCBUaGUgeCBjb21wb25lbnQgb2YgdGhlIHZlY3RvclxuICogQHByb3BlcnR5IHtudW1iZXJ9IHkgVGhlIHkgY29tcG9uZW50IG9mIHRoZSB2ZWN0b3JcbiAqL1xuXG4vKipcbiAqIENyZWF0ZSBhIG5ldyB2ZWN0b3JcbiAqIEBwYXJhbSB7bnVtYmVyfHZlY30gW3hdIFRoZSB4IGNvbXBvbmVudCBvZiB0aGUgdmVjdG9yLCBvciBhIHZlY3RvciB0byBjb3B5XG4gKiBAcGFyYW0ge251bWJlcn0gW3ldIFRoZSB5IGNvbXBvbmVudCBvZiB0aGUgdmVjdG9yXG4gKiBAcmV0dXJuIHt2ZWN9IEEgbmV3IHZlY3RvclxuICogQGV4YW1wbGUgPGNhcHRpb24+VmFyaW91cyB3YXlzIHRvIGluaXRpYWxpc2UgYSB2ZWN0b3I8L2NhcHRpb24+XG4gKiBsZXQgYSA9IHZlYygzLCAyKTsgIC8vICgzLCAyKVxuICogbGV0IGIgPSB2ZWMoNCk7ICAgICAvLyAoNCwgNClcbiAqIGxldCBjID0gdmVjKGEpOyAgICAgLy8gKDMsIDIpXG4gKiBsZXQgZCA9IHZlYygpOyAgICAgIC8vICgwLCAwKVxuICovXG5jb25zdCB2ZWMgPSAoeCwgeSkgPT4gKCF4ICYmICF5ID9cbiAgeyB4OiAwLCB5OiAwIH0gOiAodHlwZW9mIHggPT09ICdvYmplY3QnID9cbiAgICB7IHg6IHgueCB8fCAwLCB5OiB4LnkgfHwgMCB9IDogKHkgPT09IG51bGwgfHwgeSA9PT0gdW5kZWZpbmVkID9cbiAgICAgIHsgeDogeCwgeTogeCB9IDogeyB4OiB4LCB5OiB5IH0pXG4gIClcbik7XG5cbi8qKlxuICogR2V0IHRoZSBjb21wb25lbnRzIG9mIGEgdmVjdG9yIGFzIGFuIGFycmF5XG4gKiBAcGFyYW0ge3ZlY30gYSBUaGUgdmVjdG9yIHRvIGdldCBjb21wb25lbnRzIGZyb21cbiAqIEByZXR1cm4ge0FycmF5PG51bWJlcj59IFRoZSB2ZWN0b3IgY29tcG9uZW50cyBhcyBhbiBhcnJheVxuICovXG52ZWMuY29tcG9uZW50cyA9IGEgPT4gW2EueCwgYS55XTtcblxuLyoqXG4gKiBSZXR1cm4gYSB1bml0IHZlY3RvciAoMSwgMClcbiAqIEByZXR1cm4ge3ZlY30gQSB1bml0IHZlY3RvciAoMSwgMClcbiAqL1xudmVjLnV4ID0gKCkgPT4gdmVjKDEsIDApO1xuXG4vKipcbiAqIFJldHVybiBhIHVuaXQgdmVjdG9yICgwLCAxKVxuICogQHJldHVybiB7dmVjfSBBIHVuaXQgdmVjdG9yICgwLCAxKVxuICovXG52ZWMudXkgPSAoKSA9PiB2ZWMoMCwgMSk7XG5cbi8qKlxuICogQWRkIHZlY3RvcnNcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcGFyYW0ge3ZlY30gYiBWZWN0b3IgYlxuICogQHJldHVybiB7dmVjfSBhICsgYlxuICovXG52ZWMuYWRkID0gKGEsIGIpID0+ICh7IHg6IGEueCArIGIueCwgeTogYS55ICsgYi55IH0pO1xuXG4vKipcbiAqIFNjYWxlIGEgdmVjdG9yXG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHBhcmFtIHtudW1iZXJ9IGIgU2NhbGFyIGJcbiAqIEByZXR1cm4ge3ZlY30gYSAqIGJcbiAqL1xudmVjLm11bCA9IChhLCBiKSA9PiAoeyB4OiBhLnggKiBiLCB5OiBhLnkgKiBiIH0pO1xuXG4vKipcbiAqIFN1YnRyYWN0IHZlY3RvcnNcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcGFyYW0ge3ZlY30gYiBWZWN0b3IgYlxuICogQHJldHVybiB7dmVjfSBhIC0gYlxuICovXG52ZWMuc3ViID0gKGEsIGIpID0+ICh7IHg6IGEueCAtIGIueCwgeTogYS55IC0gYi55IH0pO1xuXG4vKipcbiAqIEdldCB0aGUgbGVuZ3RoIG9mIGEgdmVjdG9yXG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHJldHVybiB7bnVtYmVyfSB8YXxcbiAqL1xudmVjLmxlbiA9IGEgPT4gTWF0aC5zcXJ0KGEueCAqIGEueCArIGEueSAqIGEueSk7XG5cbi8qKlxuICogR2V0IHRoZSBsZW5ndGggb2YgYSB2ZWN0b3IgdXNpbmcgdGF4aWNhYiBnZW9tZXRyeVxuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEByZXR1cm4ge251bWJlcn0gfGF8XG4gKi9cbnZlYy5tYW5oYXR0YW4gPSBhID0+IE1hdGguYWJzKGEueCkgKyBNYXRoLmFicyhhLnkpO1xuXG4vKipcbiAqIE5vcm1hbGlzZSBhIHZlY3RvclxuICogQHBhcmFtIHt2ZWN9IGEgVGhlIHZlY3RvciB0byBub3JtYWxpc2VcbiAqIEByZXR1cm4ge3ZlY30gXmFcbiAqL1xudmVjLm5vciA9IGEgPT4ge1xuICBsZXQgbGVuID0gdmVjLmxlbihhKTtcbiAgcmV0dXJuIGxlbiA/IHsgeDogYS54IC8gbGVuLCB5OiBhLnkgLyBsZW4gfSA6IHZlYygpO1xufTtcblxuLyoqXG4gKiBHZXQgYSBkb3QgcHJvZHVjdCBvZiB2ZWN0b3JzXG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHBhcmFtIHt2ZWN9IGIgVmVjdG9yIGJcbiAqIEByZXR1cm4ge251bWJlcn0gYSDiiJkgYlxuICovXG52ZWMuZG90ID0gKGEsIGIpID0+IGEueCAqIGIueCArIGEueSAqIGIueTtcblxuLyoqXG4gKiBSb3RhdGUgYSB2ZWN0b3IgYnkgciByYWRpYW5zXG4gKiBAcGFyYW0ge3ZlY30gYSBUaGUgdmVjdG9yIHRvIHJvdGF0ZVxuICogQHBhcmFtIHtudW1iZXJ9IHIgVGhlIGFuZ2xlIHRvIHJvdGF0ZSBieSwgbWVhc3VyZWQgaW4gcmFkaWFuc1xuICogQHJldHVybiB7dmVjfSBBIHJvdGF0ZWQgdmVjdG9yXG4gKi9cbnZlYy5yb3QgPSAoYSwgcikgPT4ge1xuICBsZXQgcyA9IE1hdGguc2luKHIpLFxuICAgIGMgPSBNYXRoLmNvcyhyKTtcbiAgcmV0dXJuIHsgeDogYyAqIGEueCAtIHMgKiBhLnksIHk6IHMgKiBhLnggKyBjICogYS55IH07XG59XG5cbi8qKlxuICogQ2hlY2sgaWYgdHdvIHZlY3RvcnMgYXJlIGVxdWFsXG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHBhcmFtIHt2ZWN9IGIgVmVjdG9yIGJcbiAqIEByZXR1cm4ge2Jvb2xlYW59IFRydWUgaWYgdmVjdG9ycyBhIGFuZCBiIGFyZSBlcXVhbCwgZmFsc2Ugb3RoZXJ3aXNlXG4gKi9cbnZlYy5lcSA9IChhLCBiKSA9PiBhLnggPT09IGIueCAmJiBhLnkgPT09IGIueTtcblxuLyoqXG4gKiBHZXQgdGhlIGFuZ2xlIG9mIGEgdmVjdG9yXG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHJldHVybiB7bnVtYmVyfSBUaGUgYW5nbGUgb2YgdmVjdG9yIGEgaW4gcmFkaWFuc1xuICovXG52ZWMucmFkID0gYSA9PiBNYXRoLmF0YW4yKGEueSwgYS54KTtcblxuLyoqXG4gKiBDb3B5IGEgdmVjdG9yXG4gKiBAcGFyYW0ge3ZlY30gYSBUaGUgdmVjdG9yIHRvIGNvcHlcbiAqIEByZXR1cm4ge3ZlY30gQSBjb3B5IG9mIHZlY3RvciBhXG4gKi9cbnZlYy5jcHkgPSBhID0+IHZlYyhhKTtcblxuLyoqXG4gKiBBIGZ1bmN0aW9uIHRvIGNhbGwgb24gZWFjaCBjb21wb25lbnQgb2YgYSB2ZWN0b3JcbiAqIEBjYWxsYmFjayB2ZWN0b3JNYXBDYWxsYmFja1xuICogQHBhcmFtIHtudW1iZXJ9IHZhbHVlIFRoZSBjb21wb25lbnQgdmFsdWVcbiAqIEBwYXJhbSB7J3gnIHwgJ3knfSBsYWJlbCBUaGUgY29tcG9uZW50IGxhYmVsICh4IG9yIHkpXG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSBtYXBwZWQgY29tcG9uZW50XG4gKi9cblxuLyoqXG4gKiBDYWxsIGEgZnVuY3Rpb24gb24gZWFjaCBjb21wb25lbnQgb2YgYSB2ZWN0b3IgYW5kIGJ1aWxkIGEgbmV3IHZlY3RvciBmcm9tIHRoZSByZXN1bHRzXG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHBhcmFtIHt2ZWN0b3JNYXBDYWxsYmFja30gZiBUaGUgZnVuY3Rpb24gdG8gY2FsbCBvbiBlYWNoIGNvbXBvbmVudCBvZiB0aGUgdmVjdG9yXG4gKiBAcmV0dXJuIHt2ZWN9IFZlY3RvciBhIG1hcHBlZCB0aHJvdWdoIGZcbiAqL1xudmVjLm1hcCA9IChhLCBmKSA9PiAoeyB4OiBmKGEueCwgJ3gnKSwgeTogZihhLnksICd5JykgfSk7XG5cbi8qKlxuICogQ29udmVydCBhIHZlY3RvciBpbnRvIGEgc3RyaW5nXG4gKiBAcGFyYW0ge3ZlY30gYSBUaGUgdmVjdG9yIHRvIGNvbnZlcnRcbiAqIEBwYXJhbSB7c3RyaW5nfSBbcz0nLCAnXSBUaGUgc2VwYXJhdG9yIHN0cmluZ1xuICogQHJldHVybiB7c3RyaW5nfSBBIHN0cmluZyByZXByZXNlbnRhdGlvbiBvZiB0aGUgdmVjdG9yXG4gKi9cbnZlYy5zdHIgPSAoYSwgcyA9ICcsICcpID0+IGAke2EueH0ke3N9JHthLnl9YDtcblxuLyoqXG4gKiBBIG1hdHJpeFxuICogQHR5cGVkZWYge09iamVjdH0gbWF0XG4gKiBAcHJvcGVydHkge251bWJlcn0gbSBUaGUgbnVtYmVyIG9mIHJvd3MgaW4gdGhlIG1hdHJpeFxuICogQHByb3BlcnR5IHtudW1iZXJ9IG4gVGhlIG51bWJlciBvZiBjb2x1bW5zIGluIHRoZSBtYXRyaXhcbiAqIEBwcm9wZXJ0eSB7QXJyYXk8bnVtYmVyPn0gZW50cmllcyBUaGUgbWF0cml4IHZhbHVlc1xuICovXG5cbi8qKlxuICogQ3JlYXRlIGEgbmV3IG1hdHJpeFxuICogQHBhcmFtIHtudW1iZXJ9IFttPTRdIFRoZSBudW1iZXIgb2Ygcm93c1xuICogQHBhcmFtIHtudW1iZXJ9IFtuPTRdIFRoZSBudW1iZXIgb2YgY29sdW1uc1xuICogQHBhcmFtIHtBcnJheTxudW1iZXI+fSBbZW50cmllcz1bXV0gTWF0cml4IHZhbHVlcyBpbiByZWFkaW5nIG9yZGVyXG4gKiBAcmV0dXJuIHttYXR9IEEgbmV3IG1hdHJpeFxuICovXG5jb25zdCBtYXQgPSAobSA9IDQsIG4gPSA0LCBlbnRyaWVzID0gW10pID0+ICh7XG4gIG0sIG4sXG4gIGVudHJpZXM6IGVudHJpZXMuY29uY2F0KEFycmF5KG0gKiBuKS5maWxsKDApKS5zbGljZSgwLCBtICogbilcbn0pO1xuXG4vKipcbiAqIEdldCBhbiBpZGVudGl0eSBtYXRyaXggb2Ygc2l6ZSBuXG4gKiBAcGFyYW0ge251bWJlcn0gbiBUaGUgc2l6ZSBvZiB0aGUgbWF0cml4XG4gKiBAcmV0dXJuIHttYXR9IEFuIGlkZW50aXR5IG1hdHJpeFxuICovXG5tYXQuaWRlbnRpdHkgPSBuID0+IG1hdChuLCBuLCBBcnJheShuICogbikuZmlsbCgwKS5tYXAoKHYsIGkpID0+ICsoTWF0aC5mbG9vcihpIC8gbikgPT09IGkgJSBuKSkpO1xuXG4vKipcbiAqIEdldCBhbiBlbnRyeSBmcm9tIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIHJvdyBvZmZzZXRcbiAqIEBwYXJhbSB7bnVtYmVyfSBqIFRoZSBjb2x1bW4gb2Zmc2V0XG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSB2YWx1ZSBhdCBwb3NpdGlvbiAoaSwgaikgaW4gbWF0cml4IGFcbiAqL1xubWF0LmdldCA9IChhLCBpLCBqKSA9PiBhLmVudHJpZXNbKGogLSAxKSArIChpIC0gMSkgKiBhLm5dO1xuXG4vKipcbiAqIFNldCBhbiBlbnRyeSBvZiBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSByb3cgb2Zmc2V0XG4gKiBAcGFyYW0ge251bWJlcn0gaiBUaGUgY29sdW1uIG9mZnNldFxuICogQHBhcmFtIHtudW1iZXJ9IHYgVGhlIHZhbHVlIHRvIHNldCBpbiBtYXRyaXggYVxuICovXG5tYXQuc2V0ID0gKGEsIGksIGosIHYpID0+IHsgYS5lbnRyaWVzWyhqIC0gMSkgKyAoaSAtIDEpICogYS5uXSA9IHY7IH07XG5cbi8qKlxuICogR2V0IGEgcm93IGZyb20gYSBtYXRyaXggYXMgYW4gYXJyYXlcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge251bWJlcn0gbSBUaGUgcm93IG9mZnNldFxuICogQHJldHVybiB7QXJyYXk8bnVtYmVyPn0gUm93IG0gZnJvbSBtYXRyaXggYVxuICovXG5tYXQucm93ID0gKGEsIG0pID0+IHtcbiAgY29uc3QgcyA9IChtIC0gMSkgKiBhLm47XG4gIHJldHVybiBhLmVudHJpZXMuc2xpY2UocywgcyArIGEubik7XG59O1xuXG4vKipcbiAqIEdldCBhIGNvbHVtbiBmcm9tIGEgbWF0cml4IGFzIGFuIGFycmF5XG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHtudW1iZXJ9IG4gVGhlIGNvbHVtbiBvZmZzZXRcbiAqIEByZXR1cm4ge0FycmF5PG51bWJlcj59IENvbHVtbiBuIGZyb20gbWF0cml4IGFcbiAqL1xubWF0LmNvbCA9IChhLCBuKSA9PiB0aW1lcyhpID0+IG1hdC5nZXQoYSwgKGkgKyAxKSwgbiksIGEubSk7XG5cbi8qKlxuICogQWRkIG1hdHJpY2VzXG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHttYXR9IGIgTWF0cml4IGJcbiAqIEByZXR1cm4ge21hdH0gYSArIGJcbiAqL1xubWF0LmFkZCA9IChhLCBiKSA9PiBhLm0gPT09IGIubSAmJiBhLm4gPT09IGIubiAmJiBtYXQubWFwKGEsICh2LCBpKSA9PiB2ICsgYi5lbnRyaWVzW2ldKTtcblxuLyoqXG4gKiBTdWJ0cmFjdCBtYXRyaWNlc1xuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bWF0fSBiIE1hdHJpeCBiXG4gKiBAcmV0dXJuIHttYXR9IGEgLSBiXG4gKi9cbm1hdC5zdWIgPSAoYSwgYikgPT4gYS5tID09PSBiLm0gJiYgYS5uID09PSBiLm4gJiYgbWF0Lm1hcChhLCAodiwgaSkgPT4gdiAtIGIuZW50cmllc1tpXSk7XG5cbi8qKlxuICogTXVsdGlwbHkgbWF0cmljZXNcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge21hdH0gYiBNYXRyaXggYlxuICogQHJldHVybiB7bWF0fGJvb2xlYW59IGFiIG9yIGZhbHNlIGlmIHRoZSBtYXRyaWNlcyBjYW5ub3QgYmUgbXVsdGlwbGllZFxuICovXG5tYXQubXVsID0gKGEsIGIpID0+IHtcbiAgaWYgKGEubiAhPT0gYi5tKSB7IHJldHVybiBmYWxzZTsgfVxuICBjb25zdCByZXN1bHQgPSBtYXQoYS5tLCBiLm4pO1xuICBmb3IgKGxldCBpID0gMTsgaSA8PSBhLm07IGkrKykge1xuICAgIGZvciAobGV0IGogPSAxOyBqIDw9IGIubjsgaisrKSB7XG4gICAgICBtYXQuc2V0KHJlc3VsdCwgaSwgaiwgZG90KG1hdC5yb3coYSwgaSksIG1hdC5jb2woYiwgaikpKTtcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc3VsdDtcbn07XG5cbi8qKlxuICogU2NhbGUgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge251bWJlcn0gYiBTY2FsYXIgYlxuICogQHJldHVybiB7bWF0fSBhICogYlxuICovXG5tYXQuc2NhbGUgPSAoYSwgYikgPT4gbWF0Lm1hcChhLCB2ID0+IHYgKiBiKTtcblxuLyoqXG4gKiBUcmFuc3Bvc2UgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIFRoZSBtYXRyaXggdG8gdHJhbnNwb3NlXG4gKiBAcmV0dXJuIHttYXR9IEEgdHJhbnNwb3NlZCBtYXRyaXhcbiAqL1xubWF0LnRyYW5zID0gYSA9PiBtYXQoYS5uLCBhLm0sIHRpbWVzKGkgPT4gbWF0LmNvbChhLCAoaSArIDEpKSwgYS5uKS5mbGF0KCkpO1xuXG4vKipcbiAqIEdldCB0aGUgbWlub3Igb2YgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgcm93IG9mZnNldFxuICogQHBhcmFtIHtudW1iZXJ9IGogVGhlIGNvbHVtbiBvZmZzZXRcbiAqIEByZXR1cm4ge21hdHxib29sZWFufSBUaGUgKGksIGopIG1pbm9yIG9mIG1hdHJpeCBhIG9yIGZhbHNlIGlmIHRoZSBtYXRyaXggaXMgbm90IHNxdWFyZVxuICovXG5tYXQubWlub3IgPSAoYSwgaSwgaikgPT4ge1xuICBpZiAoYS5tICE9PSBhLm4pIHsgcmV0dXJuIGZhbHNlOyB9XG4gIGNvbnN0IGVudHJpZXMgPSBbXTtcbiAgZm9yIChsZXQgaWkgPSAxOyBpaSA8PSBhLm07IGlpKyspIHtcbiAgICBpZiAoaWkgPT09IGkpIHsgY29udGludWU7IH1cbiAgICBmb3IgKGxldCBqaiA9IDE7IGpqIDw9IGEubjsgamorKykge1xuICAgICAgaWYgKGpqID09PSBqKSB7IGNvbnRpbnVlOyB9XG4gICAgICBlbnRyaWVzLnB1c2gobWF0LmdldChhLCBpaSwgamopKTtcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG1hdChhLm0gLSAxLCBhLm4gLSAxLCBlbnRyaWVzKTtcbn07XG5cbi8qKlxuICogR2V0IHRoZSBkZXRlcm1pbmFudCBvZiBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEByZXR1cm4ge251bWJlcnxib29sZWFufSB8YXwgb3IgZmFsc2UgaWYgdGhlIG1hdHJpeCBpcyBub3Qgc3F1YXJlXG4gKi9cbm1hdC5kZXQgPSBhID0+IHtcbiAgaWYgKGEubSAhPT0gYS5uKSB7IHJldHVybiBmYWxzZTsgfVxuICBpZiAoYS5tID09PSAxKSB7XG4gICAgcmV0dXJuIGEuZW50cmllc1swXTtcbiAgfVxuICBpZiAoYS5tID09PSAyKSB7XG4gICAgcmV0dXJuIGEuZW50cmllc1swXSAqIGEuZW50cmllc1szXSAtIGEuZW50cmllc1sxXSAqIGEuZW50cmllc1syXTtcbiAgfVxuICBsZXQgdG90YWwgPSAwLCBzaWduID0gMTtcbiAgZm9yIChsZXQgaiA9IDE7IGogPD0gYS5uOyBqKyspIHtcbiAgICB0b3RhbCArPSBzaWduICogYS5lbnRyaWVzW2ogLSAxXSAqIG1hdC5kZXQobWF0Lm1pbm9yKGEsIDEsIGopKTtcbiAgICBzaWduICo9IC0xO1xuICB9XG4gIHJldHVybiB0b3RhbDtcbn07XG5cbi8qKlxuICogTm9ybWFsaXNlIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBUaGUgbWF0cml4IHRvIG5vcm1hbGlzZVxuICogQHJldHVybiB7bWF0fGJvb2xlYW59IF5hIG9yIGZhbHNlIGlmIHRoZSBtYXRyaXggaXMgbm90IHNxdWFyZVxuICovXG5tYXQubm9yID0gYSA9PiB7XG4gIGlmIChhLm0gIT09IGEubikgeyByZXR1cm4gZmFsc2U7IH1cbiAgY29uc3QgZCA9IG1hdC5kZXQoYSk7XG4gIHJldHVybiBtYXQubWFwKGEsIGkgPT4gaSAqIGQpO1xufTtcblxuLyoqXG4gKiBHZXQgdGhlIGFkanVnYXRlIG9mIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBUaGUgbWF0cml4IGZyb20gd2hpY2ggdG8gZ2V0IHRoZSBhZGp1Z2F0ZVxuICogQHJldHVybiB7bWF0fSBUaGUgYWRqdWdhdGUgb2YgYVxuICovXG5tYXQuYWRqID0gYSA9PiB7XG4gIGNvbnN0IG1pbm9ycyA9IG1hdChhLm0sIGEubik7XG4gIGZvciAobGV0IGkgPSAxOyBpIDw9IGEubTsgaSsrKSB7XG4gICAgZm9yIChsZXQgaiA9IDE7IGogPD0gYS5uOyBqKyspIHtcbiAgICAgIG1hdC5zZXQobWlub3JzLCBpLCBqLCBtYXQuZGV0KG1hdC5taW5vcihhLCBpLCBqKSkpO1xuICAgIH1cbiAgfVxuICBjb25zdCBjb2ZhY3RvcnMgPSBtYXQubWFwKG1pbm9ycywgKHYsIGkpID0+IHYgKiAoaSAlIDIgPyAtMSA6IDEpKTtcbiAgcmV0dXJuIG1hdC50cmFucyhjb2ZhY3RvcnMpO1xufTtcblxuLyoqXG4gKiBHZXQgdGhlIGludmVyc2Ugb2YgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIFRoZSBtYXRyaXggdG8gaW52ZXJ0XG4gKiBAcmV0dXJuIHttYXR8Ym9vbGVhbn0gYV4tMSBvciBmYWxzZSBpZiB0aGUgbWF0cml4IGhhcyBubyBpbnZlcnNlXG4gKi9cbm1hdC5pbnYgPSBhID0+IHtcbiAgaWYgKGEubSAhPT0gYS5uKSB7IHJldHVybiBmYWxzZTsgfVxuICBjb25zdCBkID0gbWF0LmRldChhKTtcbiAgaWYgKGQgPT09IDApIHsgcmV0dXJuIGZhbHNlOyB9XG4gIHJldHVybiBtYXQuc2NhbGUobWF0LmFkaihhKSwgMSAvIGQpO1xufTtcblxuLyoqXG4gKiBDaGVjayBpZiB0d28gbWF0cmljZXMgYXJlIGVxdWFsXG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHttYXR9IGIgTWF0cml4IGJcbiAqIEByZXR1cm4ge2Jvb2xlYW59IFRydWUgaWYgbWF0cmljZXMgYSBhbmQgYiBhcmUgaWRlbnRpY2FsLCBmYWxzZSBvdGhlcndpc2VcbiAqL1xubWF0LmVxID0gKGEsIGIpID0+IGEubSA9PT0gYi5tICYmIGEubiA9PT0gYi5uICYmIG1hdC5zdHIoYSkgPT09IG1hdC5zdHIoYik7XG5cbi8qKlxuICogQ29weSBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgVGhlIG1hdHJpeCB0byBjb3B5XG4gKiBAcmV0dXJuIHttYXR9IEEgY29weSBvZiBtYXRyaXggYVxuICovXG5tYXQuY3B5ID0gYSA9PiBtYXQoYS5tLCBhLm4sIFsuLi5hLmVudHJpZXNdKTtcblxuLyoqXG4gKiBBIGZ1bmN0aW9uIHRvIGNhbGwgb24gZWFjaCBlbnRyeSBvZiBhIG1hdHJpeFxuICogQGNhbGxiYWNrIG1hdHJpeE1hcENhbGxiYWNrXG4gKiBAcGFyYW0ge251bWJlcn0gdmFsdWUgVGhlIGVudHJ5IHZhbHVlXG4gKiBAcGFyYW0ge251bWJlcn0gaW5kZXggVGhlIGVudHJ5IGluZGV4XG4gKiBAcGFyYW0ge0FycmF5PG51bWJlcj59IGVudHJpZXMgVGhlIGFycmF5IG9mIG1hdHJpeCBlbnRyaWVzXG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSBtYXBwZWQgZW50cnlcbiAqL1xuXG4vKipcbiAqIENhbGwgYSBmdW5jdGlvbiBvbiBlYWNoIGVudHJ5IG9mIGEgbWF0cml4IGFuZCBidWlsZCBhIG5ldyBtYXRyaXggZnJvbSB0aGUgcmVzdWx0c1xuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bWF0cml4TWFwQ2FsbGJhY2t9IGYgVGhlIGZ1bmN0aW9uIHRvIGNhbGwgb24gZWFjaCBlbnRyeSBvZiB0aGUgbWF0cml4XG4gKiBAcmV0dXJuIHttYXR9IE1hdHJpeCBhIG1hcHBlZCB0aHJvdWdoIGZcbiAqL1xubWF0Lm1hcCA9IChhLCBmKSA9PiBtYXQoYS5tLCBhLm4sIGEuZW50cmllcy5tYXAoZikpO1xuXG4vKipcbiAqIENvbnZlcnQgYSBtYXRyaXggaW50byBhIHN0cmluZ1xuICogQHBhcmFtIHttYXR9IGEgVGhlIG1hdHJpeCB0byBjb252ZXJ0XG4gKiBAcGFyYW0ge3N0cmluZ30gW21zPScsICddIFRoZSBzZXBhcmF0b3Igc3RyaW5nIGZvciBjb2x1bW5zXG4gKiBAcGFyYW0ge3N0cmluZ30gW25zPSdcXG4nXSBUaGUgc2VwYXJhdG9yIHN0cmluZyBmb3Igcm93c1xuICogQHJldHVybiB7c3RyaW5nfSBBIHN0cmluZyByZXByZXNlbnRhdGlvbiBvZiB0aGUgbWF0cml4XG4gKi9cbm1hdC5zdHIgPSAoYSwgbXMgPSAnLCAnLCBucyA9ICdcXG4nKSA9PiBjaHVuayhhLmVudHJpZXMsIGEubikubWFwKHIgPT4gci5qb2luKG1zKSkuam9pbihucyk7XG5cbmlmICh0eXBlb2YgbW9kdWxlICE9PSAndW5kZWZpbmVkJykge1xuICBtb2R1bGUuZXhwb3J0cyA9IHsgdmVjLCBtYXQgfTtcbn1cbiIsImNvbnN0IGRlY29kZSA9IG51bXMgPT4ge1xuICBjb25zdCBkZWNvZGVkID0gW107XG4gIGZvciAobGV0IGkgPSAwOyBpIDwgbnVtcy5sZW5ndGg7IGkgKz0gMikge1xuICAgIGNvbnN0IHJ1bl9sZW5ndGggPSBudW1zW2ldO1xuICAgIGNvbnN0IHZhbHVlID0gbnVtc1tpICsgMV07XG4gICAgZm9yIChsZXQgaWkgPSAwOyBpaSA8IHJ1bl9sZW5ndGg7IGlpKyspIHtcbiAgICAgIGRlY29kZWQucHVzaCh2YWx1ZSk7XG4gICAgfVxuICB9XG4gIHJldHVybiBkZWNvZGVkO1xufTtcblxubW9kdWxlLmV4cG9ydHMgPSBkZWNvZGU7XG4iLCIhZnVuY3Rpb24oZyxjKXt0eXBlb2YgZXhwb3J0cz09XCJvYmplY3RcIiYmdHlwZW9mIG1vZHVsZSE9XCJ1bmRlZmluZWRcIj9jKGV4cG9ydHMpOnR5cGVvZiBkZWZpbmU9PVwiZnVuY3Rpb25cIiYmZGVmaW5lLmFtZD9kZWZpbmUoW1wiZXhwb3J0c1wiXSxjKTpjKChnPWd8fHNlbGYpLmxydV9tYXA9Zy5scnVfbWFwfHx7fSl9KHRoaXMsZnVuY3Rpb24oZyl7Y29uc3QgYz1TeW1ib2woXCJuZXdlclwiKSxlPVN5bWJvbChcIm9sZGVyXCIpO2NsYXNzIG57Y29uc3RydWN0b3IoYSxiKXt0eXBlb2YgYSE9PVwibnVtYmVyXCImJihiPWEsYT0wKSx0aGlzLnNpemU9MCx0aGlzLmxpbWl0PWEsdGhpcy5vbGRlc3Q9dGhpcy5uZXdlc3Q9dm9pZCAwLHRoaXMuX2tleW1hcD1uZXcgTWFwKCksYiYmKHRoaXMuYXNzaWduKGIpLGE8MSYmKHRoaXMubGltaXQ9dGhpcy5zaXplKSl9X21hcmtFbnRyeUFzVXNlZChhKXtpZihhPT09dGhpcy5uZXdlc3QpcmV0dXJuO2FbY10mJihhPT09dGhpcy5vbGRlc3QmJih0aGlzLm9sZGVzdD1hW2NdKSxhW2NdW2VdPWFbZV0pLGFbZV0mJihhW2VdW2NdPWFbY10pLGFbY109dm9pZCAwLGFbZV09dGhpcy5uZXdlc3QsdGhpcy5uZXdlc3QmJih0aGlzLm5ld2VzdFtjXT1hKSx0aGlzLm5ld2VzdD1hfWFzc2lnbihhKXtsZXQgYixkPXRoaXMubGltaXR8fE51bWJlci5NQVhfVkFMVUU7dGhpcy5fa2V5bWFwLmNsZWFyKCk7bGV0IG09YVtTeW1ib2wuaXRlcmF0b3JdKCk7Zm9yKGxldCBoPW0ubmV4dCgpOyFoLmRvbmU7aD1tLm5leHQoKSl7bGV0IGY9bmV3IGwoaC52YWx1ZVswXSxoLnZhbHVlWzFdKTt0aGlzLl9rZXltYXAuc2V0KGYua2V5LGYpLGI/KGJbY109ZixmW2VdPWIpOnRoaXMub2xkZXN0PWYsYj1mO2lmKGQtLT09MCl0aHJvdyBuZXcgRXJyb3IoXCJvdmVyZmxvd1wiKX10aGlzLm5ld2VzdD1iLHRoaXMuc2l6ZT10aGlzLl9rZXltYXAuc2l6ZX1nZXQoYSl7dmFyIGI9dGhpcy5fa2V5bWFwLmdldChhKTtyZXR1cm4gYj8odGhpcy5fbWFya0VudHJ5QXNVc2VkKGIpLGIudmFsdWUpOnZvaWQgMH1zZXQoYSxiKXt2YXIgZD10aGlzLl9rZXltYXAuZ2V0KGEpO3JldHVybiBkPyhkLnZhbHVlPWIsdGhpcy5fbWFya0VudHJ5QXNVc2VkKGQpLHRoaXMpOih0aGlzLl9rZXltYXAuc2V0KGEsZD1uZXcgbChhLGIpKSx0aGlzLm5ld2VzdD8odGhpcy5uZXdlc3RbY109ZCxkW2VdPXRoaXMubmV3ZXN0KTp0aGlzLm9sZGVzdD1kLHRoaXMubmV3ZXN0PWQsKyt0aGlzLnNpemUsdGhpcy5zaXplPnRoaXMubGltaXQmJnRoaXMuc2hpZnQoKSx0aGlzKX1zaGlmdCgpe3ZhciBhPXRoaXMub2xkZXN0O2lmKGEpcmV0dXJuIHRoaXMub2xkZXN0W2NdPyh0aGlzLm9sZGVzdD10aGlzLm9sZGVzdFtjXSx0aGlzLm9sZGVzdFtlXT12b2lkIDApOih0aGlzLm9sZGVzdD12b2lkIDAsdGhpcy5uZXdlc3Q9dm9pZCAwKSxhW2NdPWFbZV09dm9pZCAwLHRoaXMuX2tleW1hcC5kZWxldGUoYS5rZXkpLC0tdGhpcy5zaXplLFthLmtleSxhLnZhbHVlXX1maW5kKGEpe2xldCBiPXRoaXMuX2tleW1hcC5nZXQoYSk7cmV0dXJuIGI/Yi52YWx1ZTp2b2lkIDB9aGFzKGEpe3JldHVybiB0aGlzLl9rZXltYXAuaGFzKGEpfWRlbGV0ZShhKXt2YXIgYj10aGlzLl9rZXltYXAuZ2V0KGEpO3JldHVybiBiPyh0aGlzLl9rZXltYXAuZGVsZXRlKGIua2V5KSxiW2NdJiZiW2VdPyhiW2VdW2NdPWJbY10sYltjXVtlXT1iW2VdKTpiW2NdPyhiW2NdW2VdPXZvaWQgMCx0aGlzLm9sZGVzdD1iW2NdKTpiW2VdPyhiW2VdW2NdPXZvaWQgMCx0aGlzLm5ld2VzdD1iW2VdKTp0aGlzLm9sZGVzdD10aGlzLm5ld2VzdD12b2lkIDAsdGhpcy5zaXplLS0sYi52YWx1ZSk6dm9pZCAwfWNsZWFyKCl7dGhpcy5vbGRlc3Q9dGhpcy5uZXdlc3Q9dm9pZCAwLHRoaXMuc2l6ZT0wLHRoaXMuX2tleW1hcC5jbGVhcigpfWtleXMoKXtyZXR1cm4gbmV3IGoodGhpcy5vbGRlc3QpfXZhbHVlcygpe3JldHVybiBuZXcgayh0aGlzLm9sZGVzdCl9ZW50cmllcygpe3JldHVybiB0aGlzfVtTeW1ib2wuaXRlcmF0b3JdKCl7cmV0dXJuIG5ldyBpKHRoaXMub2xkZXN0KX1mb3JFYWNoKGEsYil7dHlwZW9mIGIhPT1cIm9iamVjdFwiJiYoYj10aGlzKTtsZXQgZD10aGlzLm9sZGVzdDtmb3IoO2Q7KWEuY2FsbChiLGQudmFsdWUsZC5rZXksdGhpcyksZD1kW2NdfXRvSlNPTigpe2Zvcih2YXIgYT1uZXcgQXJyYXkodGhpcy5zaXplKSxiPTAsZD10aGlzLm9sZGVzdDtkOylhW2IrK109e2tleTpkLmtleSx2YWx1ZTpkLnZhbHVlfSxkPWRbY107cmV0dXJuIGF9dG9TdHJpbmcoKXtmb3IodmFyIGE9XCJcIixiPXRoaXMub2xkZXN0O2I7KWErPVN0cmluZyhiLmtleSkrXCI6XCIrYi52YWx1ZSxiPWJbY10sYiYmKGErPVwiIDwgXCIpO3JldHVybiBhfX1nLkxSVU1hcD1uO2Z1bmN0aW9uIGwoYSxiKXt0aGlzLmtleT1hLHRoaXMudmFsdWU9Yix0aGlzW2NdPXZvaWQgMCx0aGlzW2VdPXZvaWQgMH1mdW5jdGlvbiBpKGEpe3RoaXMuZW50cnk9YX1pLnByb3RvdHlwZVtTeW1ib2wuaXRlcmF0b3JdPWZ1bmN0aW9uKCl7cmV0dXJuIHRoaXN9LGkucHJvdG90eXBlLm5leHQ9ZnVuY3Rpb24oKXtsZXQgYT10aGlzLmVudHJ5O3JldHVybiBhPyh0aGlzLmVudHJ5PWFbY10se2RvbmU6ITEsdmFsdWU6W2Eua2V5LGEudmFsdWVdfSk6e2RvbmU6ITAsdmFsdWU6dm9pZCAwfX07ZnVuY3Rpb24gaihhKXt0aGlzLmVudHJ5PWF9ai5wcm90b3R5cGVbU3ltYm9sLml0ZXJhdG9yXT1mdW5jdGlvbigpe3JldHVybiB0aGlzfSxqLnByb3RvdHlwZS5uZXh0PWZ1bmN0aW9uKCl7bGV0IGE9dGhpcy5lbnRyeTtyZXR1cm4gYT8odGhpcy5lbnRyeT1hW2NdLHtkb25lOiExLHZhbHVlOmEua2V5fSk6e2RvbmU6ITAsdmFsdWU6dm9pZCAwfX07ZnVuY3Rpb24gayhhKXt0aGlzLmVudHJ5PWF9ay5wcm90b3R5cGVbU3ltYm9sLml0ZXJhdG9yXT1mdW5jdGlvbigpe3JldHVybiB0aGlzfSxrLnByb3RvdHlwZS5uZXh0PWZ1bmN0aW9uKCl7bGV0IGE9dGhpcy5lbnRyeTtyZXR1cm4gYT8odGhpcy5lbnRyeT1hW2NdLHtkb25lOiExLHZhbHVlOmEudmFsdWV9KTp7ZG9uZTohMCx2YWx1ZTp2b2lkIDB9fX0pO1xuLy8jIHNvdXJjZU1hcHBpbmdVUkw9bHJ1LmpzLm1hcFxuIiwiXCJ1c2Ugc3RyaWN0XCI7XG5PYmplY3QuZGVmaW5lUHJvcGVydHkoZXhwb3J0cywgXCJfX2VzTW9kdWxlXCIsIHsgdmFsdWU6IHRydWUgfSk7XG5leHBvcnRzLmJpdG1hcFRvUmVjdGFuZ2xlcyA9IHZvaWQgMDtcbmNvbnN0IHZlY18xID0gcmVxdWlyZShcIkBiYXNlbWVudHVuaXZlcnNlL3ZlY1wiKTtcbmZ1bmN0aW9uIGJpdG1hcFRvUmVjdGFuZ2xlcyhiaXRtYXApIHtcbiAgICBjb25zdCByZWN0YW5nbGVzID0gW107XG4gICAgLy8gU3RlcCAxIC0gY3JlYXRlIDEtdW5pdCB0YWxsIHJlY3RhbmdsZXMgZm9yIGVhY2ggcm93XG4gICAgZm9yIChjb25zdCBbeSwgcm93XSBvZiBiaXRtYXAuZW50cmllcygpKSB7XG4gICAgICAgIGxldCBjdXJyZW50UmVjdGFuZ2xlID0gbnVsbDtcbiAgICAgICAgZm9yIChsZXQgeCA9IDA7IHggPCByb3cubGVuZ3RoOyB4KyspIHtcbiAgICAgICAgICAgIGlmIChyb3dbeF0pIHtcbiAgICAgICAgICAgICAgICBpZiAoIWN1cnJlbnRSZWN0YW5nbGUpIHtcbiAgICAgICAgICAgICAgICAgICAgY3VycmVudFJlY3RhbmdsZSA9IHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHBvc2l0aW9uOiAoMCwgdmVjXzEudmVjKSh4LCB5KSxcbiAgICAgICAgICAgICAgICAgICAgICAgIHNpemU6ICgwLCB2ZWNfMS52ZWMpKDEsIDEpLFxuICAgICAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgY3VycmVudFJlY3RhbmdsZS5zaXplLngrKztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgICAgICBpZiAoY3VycmVudFJlY3RhbmdsZSkge1xuICAgICAgICAgICAgICAgICAgICByZWN0YW5nbGVzLnB1c2goY3VycmVudFJlY3RhbmdsZSk7XG4gICAgICAgICAgICAgICAgICAgIGN1cnJlbnRSZWN0YW5nbGUgPSBudWxsO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cbiAgICAvLyBTdGVwIDIgLSBleHRlbmQgZWFjaCByZWN0YW5nbGUgZG93bndhcmRzIGlmIHBvc3NpYmxlXG4gICAgbGV0IHBhaXI7XG4gICAgd2hpbGUgKHBhaXIgPSBmaW5kUmVjdGFuZ2xlVG9FeHRlbmQocmVjdGFuZ2xlcykpIHtcbiAgICAgICAgY29uc3QgW2EsIGJdID0gcGFpcjtcbiAgICAgICAgcmVjdGFuZ2xlcy5zcGxpY2UoaW5kZXhPZihiLCByZWN0YW5nbGVzKSwgMSwgLi4uY2hvcFJlY3RhbmdsZShiLCBhKSk7XG4gICAgICAgIGEuc2l6ZS55ICs9IGIuc2l6ZS55O1xuICAgIH1cbiAgICByZXR1cm4gcmVjdGFuZ2xlcztcbn1cbmV4cG9ydHMuYml0bWFwVG9SZWN0YW5nbGVzID0gYml0bWFwVG9SZWN0YW5nbGVzO1xuLyoqXG4gKiBHZXQgdGhlIGluZGV4IG9mIHJlY3RhbmdsZSBhIGluIGEgbGlzdCBvZiByZWN0YW5nbGVzXG4gKi9cbmZ1bmN0aW9uIGluZGV4T2YoYSwgcmVjdGFuZ2xlcykge1xuICAgIHJldHVybiByZWN0YW5nbGVzLmZpbmRJbmRleChiID0+IHZlY18xLnZlYy5lcShhLnBvc2l0aW9uLCBiLnBvc2l0aW9uKSAmJiB2ZWNfMS52ZWMuZXEoYS5zaXplLCBiLnNpemUpKTtcbn1cbi8qKlxuICogRmluZCBhIHBhaXIgb2YgcmVjdGFuZ2xlcyB3aGVyZSB0aGUgZmlyc3Qgb25lIGNhbiBiZSBleHRlbmRlZCBpbnRvIHRoZVxuICogc2Vjb25kIG9uZVxuICpcbiAqIElmIG5vIHN1Y2ggcGFpciBleGlzdHMsIHJldHVybiBudWxsXG4gKi9cbmZ1bmN0aW9uIGZpbmRSZWN0YW5nbGVUb0V4dGVuZChyZWN0YW5nbGVzKSB7XG4gICAgZm9yIChjb25zdCBhIG9mIHJlY3RhbmdsZXMpIHtcbiAgICAgICAgY29uc3QgYiA9IGZpbmRSZWN0YW5nbGVUb0V4dGVuZEludG8oYSwgcmVjdGFuZ2xlcyk7XG4gICAgICAgIGlmIChiKSB7XG4gICAgICAgICAgICByZXR1cm4gW2EsIGJdO1xuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiBudWxsO1xufVxuLyoqXG4gKiBGaW5kIGEgcmVjdGFuZ2xlIHdoaWNoIHJlY3RhbmdsZSBhIGNhbiBleHRlbmQgaW50bywgb3IgbnVsbCBpZiBub25lIGNhbiBiZVxuICogZm91bmRcbiAqXG4gKiBBIHJlY3RhbmdsZSBjYW4gZXh0ZW5kIGludG8gYW5vdGhlciBvbmUgaWYgdGhlIG90aGVyIG9uZSBpcyBleGFjdGx5IGJlbG93XG4gKiBhbmQgdGhlaXIgeC1heGlzIHByb2plY3Rpb25zIG92ZXJsYXBcbiAqL1xuZnVuY3Rpb24gZmluZFJlY3RhbmdsZVRvRXh0ZW5kSW50byhhLCByZWN0YW5nbGVzKSB7XG4gICAgdmFyIF9hO1xuICAgIHJldHVybiAoX2EgPSByZWN0YW5nbGVzLmZpbmQob3RoZXIgPT4gKFxuICAgIC8vIFRoZSBvdGhlciByZWN0YW5nbGUgaXMgZXhhY3RseSBiZWxvdyB0aGUgY3VycmVudCBvbmVcbiAgICBvdGhlci5wb3NpdGlvbi55ID09PSBhLnBvc2l0aW9uLnkgKyBhLnNpemUueSAmJlxuICAgICAgICAvLyBUaGUgb3RoZXIgcmVjdGFuZ2xlIHN0YXJ0cyBiZWZvcmUgKG9yIGF0KSB0aGUgc3RhcnQgb2YgdGhlIGN1cnJlbnQgb25lXG4gICAgICAgIG90aGVyLnBvc2l0aW9uLnggPD0gYS5wb3NpdGlvbi54ICYmXG4gICAgICAgIC8vIFRoZSBvdGhlciByZWN0YW5nbGUgZW5kcyBhZnRlciAob3IgYXQpIHRoZSBlbmQgb2YgdGhlIGN1cnJlbnQgb25lXG4gICAgICAgIG90aGVyLnBvc2l0aW9uLnggKyBvdGhlci5zaXplLnggPj0gYS5wb3NpdGlvbi54ICsgYS5zaXplLngpKSkgIT09IG51bGwgJiYgX2EgIT09IHZvaWQgMCA/IF9hIDogbnVsbDtcbn1cbi8qKlxuICogU3VidHJhY3QgcmVjdGFuZ2xlIGIgZnJvbSByZWN0YW5nbGUgYSwgaWdub3JpbmcgaGVpZ2h0IChpLmUuIG9ubHkgaW4gdGhlXG4gKiB4LWF4aXMpIGFuZCByZXR1cm4gMCwgMSBvciAyIHJlc3VsdGluZyByZWN0YW5nbGVzXG4gKi9cbmZ1bmN0aW9uIGNob3BSZWN0YW5nbGUoYSwgYikge1xuICAgIGNvbnN0IHJlc3VsdCA9IFtdO1xuICAgIGlmIChiLnBvc2l0aW9uLnggPiBhLnBvc2l0aW9uLngpIHtcbiAgICAgICAgcmVzdWx0LnB1c2goe1xuICAgICAgICAgICAgcG9zaXRpb246ICgwLCB2ZWNfMS52ZWMpKGEucG9zaXRpb24ueCwgYS5wb3NpdGlvbi55KSxcbiAgICAgICAgICAgIHNpemU6ICgwLCB2ZWNfMS52ZWMpKGIucG9zaXRpb24ueCAtIGEucG9zaXRpb24ueCwgYS5zaXplLnkpLFxuICAgICAgICB9KTtcbiAgICB9XG4gICAgaWYgKGIucG9zaXRpb24ueCArIGIuc2l6ZS54IDwgYS5wb3NpdGlvbi54ICsgYS5zaXplLngpIHtcbiAgICAgICAgcmVzdWx0LnB1c2goe1xuICAgICAgICAgICAgcG9zaXRpb246ICgwLCB2ZWNfMS52ZWMpKGIucG9zaXRpb24ueCArIGIuc2l6ZS54LCBhLnBvc2l0aW9uLnkpLFxuICAgICAgICAgICAgc2l6ZTogKDAsIHZlY18xLnZlYykoKGEucG9zaXRpb24ueCArIGEuc2l6ZS54KSAtIChiLnBvc2l0aW9uLnggKyBiLnNpemUueCksIGEuc2l6ZS55KSxcbiAgICAgICAgfSk7XG4gICAgfVxuICAgIHJldHVybiByZXN1bHQ7XG59XG4vLyMgc291cmNlTWFwcGluZ1VSTD1kYXRhOmFwcGxpY2F0aW9uL2pzb247YmFzZTY0LGV5SjJaWEp6YVc5dUlqb3pMQ0ptYVd4bElqb2lZbWwwYldGd0xXUmxZMjl0Y0c5elpTNXFjeUlzSW5OdmRYSmpaVkp2YjNRaU9pSWlMQ0p6YjNWeVkyVnpJanBiSWk0dUwySnBkRzFoY0Mxa1pXTnZiWEJ2YzJVdWRITWlYU3dpYm1GdFpYTWlPbHRkTENKdFlYQndhVzVuY3lJNklqczdPMEZCUVVFc0swTkJRVFJETzBGQlR6VkRMRk5CUVdkQ0xHdENRVUZyUWl4RFFVRkRMRTFCUVcxQ08wbEJRM0JFTEUxQlFVMHNWVUZCVlN4SFFVRm5RaXhGUVVGRkxFTkJRVU03U1VGRmJrTXNjMFJCUVhORU8wbEJRM1JFTEV0QlFVc3NUVUZCVFN4RFFVRkRMRU5CUVVNc1JVRkJSU3hIUVVGSExFTkJRVU1zU1VGQlNTeE5RVUZOTEVOQlFVTXNUMEZCVHl4RlFVRkZMRVZCUVVVN1VVRkRka01zU1VGQlNTeG5Ra0ZCWjBJc1IwRkJjVUlzU1VGQlNTeERRVUZETzFGQlJUbERMRXRCUVVzc1NVRkJTU3hEUVVGRExFZEJRVWNzUTBGQlF5eEZRVUZGTEVOQlFVTXNSMEZCUnl4SFFVRkhMRU5CUVVNc1RVRkJUU3hGUVVGRkxFTkJRVU1zUlVGQlJTeEZRVUZGTzFsQlEyNURMRWxCUVVrc1IwRkJSeXhEUVVGRExFTkJRVU1zUTBGQlF5eEZRVUZGTzJkQ1FVTldMRWxCUVVrc1EwRkJReXhuUWtGQlowSXNSVUZCUlR0dlFrRkRja0lzWjBKQlFXZENMRWRCUVVjN2QwSkJRMnBDTEZGQlFWRXNSVUZCUlN4SlFVRkJMRk5CUVVjc1JVRkJReXhEUVVGRExFVkJRVVVzUTBGQlF5eERRVUZETzNkQ1FVTnVRaXhKUVVGSkxFVkJRVVVzU1VGQlFTeFRRVUZITEVWQlFVTXNRMEZCUXl4RlFVRkZMRU5CUVVNc1EwRkJRenR4UWtGRGFFSXNRMEZCUXp0cFFrRkRTRHR4UWtGQlRUdHZRa0ZEVEN4blFrRkJaMElzUTBGQlF5eEpRVUZKTEVOQlFVTXNRMEZCUXl4RlFVRkZMRU5CUVVNN2FVSkJRek5DTzJGQlEwWTdhVUpCUVUwN1owSkJRMHdzU1VGQlNTeG5Ra0ZCWjBJc1JVRkJSVHR2UWtGRGNFSXNWVUZCVlN4RFFVRkRMRWxCUVVrc1EwRkJReXhuUWtGQlowSXNRMEZCUXl4RFFVRkRPMjlDUVVOc1F5eG5Ra0ZCWjBJc1IwRkJSeXhKUVVGSkxFTkJRVU03YVVKQlEzcENPMkZCUTBZN1UwRkRSanRMUVVOR08wbEJSVVFzZFVSQlFYVkVPMGxCUTNaRUxFbEJRVWtzU1VGQmJVTXNRMEZCUXp0SlFVTjRReXhQUVVGUExFbEJRVWtzUjBGQlJ5eHhRa0ZCY1VJc1EwRkJReXhWUVVGVkxFTkJRVU1zUlVGQlJUdFJRVU12UXl4TlFVRk5MRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zUTBGQlF5eEhRVUZITEVsQlFVa3NRMEZCUXp0UlFVVndRaXhWUVVGVkxFTkJRVU1zVFVGQlRTeERRVUZETEU5QlFVOHNRMEZCUXl4RFFVRkRMRVZCUVVVc1ZVRkJWU3hEUVVGRExFVkJRVVVzUTBGQlF5eEZRVUZGTEVkQlFVY3NZVUZCWVN4RFFVRkRMRU5CUVVNc1JVRkJSU3hEUVVGRExFTkJRVU1zUTBGQlF5eERRVUZETzFGQlJYSkZMRU5CUVVNc1EwRkJReXhKUVVGSkxFTkJRVU1zUTBGQlF5eEpRVUZKTEVOQlFVTXNRMEZCUXl4SlFVRkpMRU5CUVVNc1EwRkJReXhEUVVGRE8wdEJRM1JDTzBsQlJVUXNUMEZCVHl4VlFVRlZMRU5CUVVNN1FVRkRjRUlzUTBGQlF6dEJRWEpEUkN4blJFRnhRME03UVVGRlJEczdSMEZGUnp0QlFVTklMRk5CUVZNc1QwRkJUeXhEUVVGRExFTkJRVmtzUlVGQlJTeFZRVUYxUWp0SlFVTndSQ3hQUVVGUExGVkJRVlVzUTBGQlF5eFRRVUZUTEVOQlEzcENMRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zVTBGQlJ5eERRVUZETEVWQlFVVXNRMEZCUXl4RFFVRkRMRU5CUVVNc1VVRkJVU3hGUVVGRkxFTkJRVU1zUTBGQlF5eFJRVUZSTEVOQlFVTXNTVUZCU1N4VFFVRkhMRU5CUVVNc1JVRkJSU3hEUVVGRExFTkJRVU1zUTBGQlF5eEpRVUZKTEVWQlFVVXNRMEZCUXl4RFFVRkRMRWxCUVVrc1EwRkJReXhEUVVNNVJDeERRVUZETzBGQlEwb3NRMEZCUXp0QlFVVkVPenM3T3p0SFFVdEhPMEZCUTBnc1UwRkJVeXh4UWtGQmNVSXNRMEZETlVJc1ZVRkJkVUk3U1VGRmRrSXNTMEZCU3l4TlFVRk5MRU5CUVVNc1NVRkJTU3hWUVVGVkxFVkJRVVU3VVVGRE1VSXNUVUZCVFN4RFFVRkRMRWRCUVVjc2VVSkJRWGxDTEVOQlFVTXNRMEZCUXl4RlFVRkZMRlZCUVZVc1EwRkJReXhEUVVGRE8xRkJRMjVFTEVsQlFVa3NRMEZCUXl4RlFVRkZPMWxCUTB3c1QwRkJUeXhEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVOQlFVTXNRMEZCUXp0VFFVTm1PMHRCUTBZN1NVRkZSQ3hQUVVGUExFbEJRVWtzUTBGQlF6dEJRVU5rTEVOQlFVTTdRVUZGUkRzN096czdPMGRCVFVjN1FVRkRTQ3hUUVVGVExIbENRVUY1UWl4RFFVTm9ReXhEUVVGWkxFVkJRMW9zVlVGQmRVSTdPMGxCUlhaQ0xFOUJRVThzVFVGQlFTeFZRVUZWTEVOQlFVTXNTVUZCU1N4RFFVTndRaXhMUVVGTExFTkJRVU1zUlVGQlJTeERRVUZETzBsQlExQXNkVVJCUVhWRU8wbEJRM1pFTEV0QlFVc3NRMEZCUXl4UlFVRlJMRU5CUVVNc1EwRkJReXhMUVVGTExFTkJRVU1zUTBGQlF5eFJRVUZSTEVOQlFVTXNRMEZCUXl4SFFVRkhMRU5CUVVNc1EwRkJReXhKUVVGSkxFTkJRVU1zUTBGQlF6dFJRVVUxUXl4NVJVRkJlVVU3VVVGRGVrVXNTMEZCU3l4RFFVRkRMRkZCUVZFc1EwRkJReXhEUVVGRExFbEJRVWtzUTBGQlF5eERRVUZETEZGQlFWRXNRMEZCUXl4RFFVRkRPMUZCUldoRExHOUZRVUZ2UlR0UlFVTndSU3hMUVVGTExFTkJRVU1zVVVGQlVTeERRVUZETEVOQlFVTXNSMEZCUnl4TFFVRkxMRU5CUVVNc1NVRkJTU3hEUVVGRExFTkJRVU1zU1VGQlNTeERRVUZETEVOQlFVTXNVVUZCVVN4RFFVRkRMRU5CUVVNc1IwRkJSeXhEUVVGRExFTkJRVU1zU1VGQlNTeERRVUZETEVOQlFVTXNRMEZETTBRc1EwRkRSaXh0UTBGQlNTeEpRVUZKTEVOQlFVTTdRVUZEV2l4RFFVRkRPMEZCUlVRN096dEhRVWRITzBGQlEwZ3NVMEZCVXl4aFFVRmhMRU5CUTNCQ0xFTkJRVmtzUlVGRFdpeERRVUZaTzBsQlJWb3NUVUZCVFN4TlFVRk5MRWRCUVdkQ0xFVkJRVVVzUTBGQlF6dEpRVU12UWl4SlFVRkpMRU5CUVVNc1EwRkJReXhSUVVGUkxFTkJRVU1zUTBGQlF5eEhRVUZITEVOQlFVTXNRMEZCUXl4UlFVRlJMRU5CUVVNc1EwRkJReXhGUVVGRk8xRkJReTlDTEUxQlFVMHNRMEZCUXl4SlFVRkpMRU5CUVVNN1dVRkRWaXhSUVVGUkxFVkJRVVVzU1VGQlFTeFRRVUZITEVWQlFVTXNRMEZCUXl4RFFVRkRMRkZCUVZFc1EwRkJReXhEUVVGRExFVkJRVVVzUTBGQlF5eERRVUZETEZGQlFWRXNRMEZCUXl4RFFVRkRMRU5CUVVNN1dVRkRla01zU1VGQlNTeEZRVUZGTEVsQlFVRXNVMEZCUnl4RlFVRkRMRU5CUVVNc1EwRkJReXhSUVVGUkxFTkJRVU1zUTBGQlF5eEhRVUZITEVOQlFVTXNRMEZCUXl4UlFVRlJMRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zUTBGQlF5eEpRVUZKTEVOQlFVTXNRMEZCUXl4RFFVRkRPMU5CUTJwRUxFTkJRVU1zUTBGQlF6dExRVU5LTzBsQlEwUXNTVUZCU1N4RFFVRkRMRU5CUVVNc1VVRkJVU3hEUVVGRExFTkJRVU1zUjBGQlJ5eERRVUZETEVOQlFVTXNTVUZCU1N4RFFVRkRMRU5CUVVNc1IwRkJSeXhEUVVGRExFTkJRVU1zVVVGQlVTeERRVUZETEVOQlFVTXNSMEZCUnl4RFFVRkRMRU5CUVVNc1NVRkJTU3hEUVVGRExFTkJRVU1zUlVGQlJUdFJRVU55UkN4TlFVRk5MRU5CUVVNc1NVRkJTU3hEUVVGRE8xbEJRMVlzVVVGQlVTeEZRVUZGTEVsQlFVRXNVMEZCUnl4RlFVRkRMRU5CUVVNc1EwRkJReXhSUVVGUkxFTkJRVU1zUTBGQlF5eEhRVUZITEVOQlFVTXNRMEZCUXl4SlFVRkpMRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zUTBGQlF5eFJRVUZSTEVOQlFVTXNRMEZCUXl4RFFVRkRPMWxCUTNCRUxFbEJRVWtzUlVGQlJTeEpRVUZCTEZOQlFVY3NSVUZEVUN4RFFVRkRMRU5CUVVNc1EwRkJReXhSUVVGUkxFTkJRVU1zUTBGQlF5eEhRVUZITEVOQlFVTXNRMEZCUXl4SlFVRkpMRU5CUVVNc1EwRkJReXhEUVVGRExFZEJRVWNzUTBGQlF5eERRVUZETEVOQlFVTXNVVUZCVVN4RFFVRkRMRU5CUVVNc1IwRkJSeXhEUVVGRExFTkJRVU1zU1VGQlNTeERRVUZETEVOQlFVTXNRMEZCUXl4RlFVTnlSQ3hEUVVGRExFTkJRVU1zU1VGQlNTeERRVUZETEVOQlFVTXNRMEZEVkR0VFFVTkdMRU5CUVVNc1EwRkJRenRMUVVOS08wbEJSVVFzVDBGQlR5eE5RVUZOTEVOQlFVTTdRVUZEYUVJc1EwRkJReUo5IiwiXCJ1c2Ugc3RyaWN0XCI7XG52YXIgX19pbXBvcnREZWZhdWx0ID0gKHRoaXMgJiYgdGhpcy5fX2ltcG9ydERlZmF1bHQpIHx8IGZ1bmN0aW9uIChtb2QpIHtcbiAgICByZXR1cm4gKG1vZCAmJiBtb2QuX19lc01vZHVsZSkgPyBtb2QgOiB7IFwiZGVmYXVsdFwiOiBtb2QgfTtcbn07XG5PYmplY3QuZGVmaW5lUHJvcGVydHkoZXhwb3J0cywgXCJfX2VzTW9kdWxlXCIsIHsgdmFsdWU6IHRydWUgfSk7XG5leHBvcnRzLnRpbGVNYXBPcHRpb25zQ29udGVudFByb2Nlc3NvciA9IGV4cG9ydHMuVGlsZU1hcCA9IGV4cG9ydHMuVGlsZUFsaWdubWVudCA9IHZvaWQgMDtcbmNvbnN0IGxydV9tYXBfMSA9IHJlcXVpcmUoXCJscnVfbWFwXCIpO1xuY29uc3QgZGVjb2RlXzEgPSBfX2ltcG9ydERlZmF1bHQocmVxdWlyZShcImZhc3QtcmxlL2RlY29kZVwiKSk7XG5jb25zdCB2ZWNfMSA9IHJlcXVpcmUoXCJAYmFzZW1lbnR1bml2ZXJzZS92ZWNcIik7XG5jb25zdCB1dGlsc18xID0gcmVxdWlyZShcIkBiYXNlbWVudHVuaXZlcnNlL3V0aWxzXCIpO1xuY29uc3QgYml0bWFwX2RlY29tcG9zZV8xID0gcmVxdWlyZShcIi4vYml0bWFwLWRlY29tcG9zZVwiKTtcbnZhciBUaWxlQWxpZ25tZW50O1xuKGZ1bmN0aW9uIChUaWxlQWxpZ25tZW50KSB7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiVG9wTGVmdFwiXSA9IDBdID0gXCJUb3BMZWZ0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiVG9wXCJdID0gMV0gPSBcIlRvcFwiO1xuICAgIFRpbGVBbGlnbm1lbnRbVGlsZUFsaWdubWVudFtcIlRvcFJpZ2h0XCJdID0gMl0gPSBcIlRvcFJpZ2h0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiTGVmdFwiXSA9IDNdID0gXCJMZWZ0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiQ2VudGVyXCJdID0gNF0gPSBcIkNlbnRlclwiO1xuICAgIFRpbGVBbGlnbm1lbnRbVGlsZUFsaWdubWVudFtcIlJpZ2h0XCJdID0gNV0gPSBcIlJpZ2h0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiQm90dG9tTGVmdFwiXSA9IDZdID0gXCJCb3R0b21MZWZ0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiQm90dG9tXCJdID0gN10gPSBcIkJvdHRvbVwiO1xuICAgIFRpbGVBbGlnbm1lbnRbVGlsZUFsaWdubWVudFtcIkJvdHRvbVJpZ2h0XCJdID0gOF0gPSBcIkJvdHRvbVJpZ2h0XCI7XG59KShUaWxlQWxpZ25tZW50ID0gZXhwb3J0cy5UaWxlQWxpZ25tZW50IHx8IChleHBvcnRzLlRpbGVBbGlnbm1lbnQgPSB7fSkpO1xuZnVuY3Rpb24gcG9pbnRJblJlY3RhbmdsZShwb2ludCwgdG9wTGVmdCwgYm90dG9tUmlnaHQpIHtcbiAgICByZXR1cm4gKHBvaW50LnggPj0gdG9wTGVmdC54ICYmXG4gICAgICAgIHBvaW50LnkgPj0gdG9wTGVmdC55ICYmXG4gICAgICAgIHBvaW50LnggPCBib3R0b21SaWdodC54ICYmXG4gICAgICAgIHBvaW50LnkgPCBib3R0b21SaWdodC55KTtcbn1cbmNsYXNzIFRpbGVNYXAge1xuICAgIGNvbnN0cnVjdG9yKG9wdGlvbnMpIHtcbiAgICAgICAgY29uc3QgYWN0dWFsT3B0aW9ucyA9IE9iamVjdC5hc3NpZ24oe30sIFRpbGVNYXAuREVGQVVMVF9PUFRJT05TLCBvcHRpb25zICE9PSBudWxsICYmIG9wdGlvbnMgIT09IHZvaWQgMCA/IG9wdGlvbnMgOiB7fSk7XG4gICAgICAgIGlmICghYWN0dWFsT3B0aW9ucy5kZWJ1ZyB8fCBhY3R1YWxPcHRpb25zLmRlYnVnID09PSB0cnVlKSB7XG4gICAgICAgICAgICBhY3R1YWxPcHRpb25zLmRlYnVnID0ge1xuICAgICAgICAgICAgICAgIHNob3dPcmlnaW46ICEhYWN0dWFsT3B0aW9ucy5kZWJ1ZyxcbiAgICAgICAgICAgICAgICBzaG93Q2h1bmtCb3JkZXJzOiAhIWFjdHVhbE9wdGlvbnMuZGVidWcsXG4gICAgICAgICAgICAgICAgc2hvd0NodW5rTGFiZWxzOiAhIWFjdHVhbE9wdGlvbnMuZGVidWcsXG4gICAgICAgICAgICAgICAgc2hvd1RpbGVCb3JkZXJzOiAhIWFjdHVhbE9wdGlvbnMuZGVidWcsXG4gICAgICAgICAgICB9O1xuICAgICAgICB9XG4gICAgICAgIHRoaXMub3B0aW9ucyA9IGFjdHVhbE9wdGlvbnM7XG4gICAgICAgIHRoaXMuY2h1bmtCdWZmZXIgPSBuZXcgbHJ1X21hcF8xLkxSVU1hcCh0aGlzLm9wdGlvbnMuY2h1bmtCdWZmZXJNYXhTaXplKTtcbiAgICB9XG4gICAgLyoqXG4gICAgICogR2V0IGEgKHJvdWdobHkgbWluaW1hbCkgc2V0IG9mIHJlY3RhbmdsZXMgd2hpY2ggY292ZXIgdGhlIHRpbGVzIGluIGFcbiAgICAgKiBnaXZlbiBsYXllclxuICAgICAqXG4gICAgICogQHBhcmFtIGxheWVyTmFtZSBUaGUgbmFtZSBvZiB0aGUgbGF5ZXIgdG8gZ2V0IHJlY3RhbmdsZXMgZm9yXG4gICAgICogQHBhcmFtIGZpZWxkTmFtZSBXZSB3aWxsIGNoZWNrIHRoZSB0cnV0aHluZXNzIG9mIHRoaXMgZmllbGQgaW4gdGhlXG4gICAgICogdGlsZSBkZWZpbml0aW9uXG4gICAgICogQHBhcmFtIHRpbGVCb3VuZHMgT3B0aW9uYWwgYm91bmRzIHRvIGNoZWNrIHdpdGhpbiwgcmVsYXRpdmUgdG8gYm91bmRzXG4gICAgICogZGVmaW5lZCBpbiBvcHRpb25zIGlmIGFueSBleGlzdCwgb3RoZXJ3aXNlIHJlbGF0aXZlIHRvICgwLCAwKVxuICAgICAqL1xuICAgIGdldExheWVyUmVjdGFuZ2xlcyhsYXllck5hbWUsIGZpZWxkTmFtZSwgdGlsZUJvdW5kcykge1xuICAgICAgICB2YXIgX2EsIF9iLCBfYywgX2QsIF9lLCBfZiwgX2csIF9oLCBfajtcbiAgICAgICAgY29uc3QgbGF5ZXIgPSB0aGlzLm9wdGlvbnMubGF5ZXJzLmZpbmQoKGwpID0+IGwubmFtZSA9PT0gbGF5ZXJOYW1lKTtcbiAgICAgICAgaWYgKCFsYXllcikge1xuICAgICAgICAgICAgcmV0dXJuIFtdO1xuICAgICAgICB9XG4gICAgICAgIGNvbnN0IHRvcExlZnQgPSAoX2EgPSB0aWxlQm91bmRzID09PSBudWxsIHx8IHRpbGVCb3VuZHMgPT09IHZvaWQgMCA/IHZvaWQgMCA6IHRpbGVCb3VuZHMudG9wTGVmdCkgIT09IG51bGwgJiYgX2EgIT09IHZvaWQgMCA/IF9hIDogKDAsIHZlY18xLnZlYykoMCk7XG4gICAgICAgIGNvbnN0IGJvdHRvbVJpZ2h0ID0gKF9iID0gdGlsZUJvdW5kcyA9PT0gbnVsbCB8fCB0aWxlQm91bmRzID09PSB2b2lkIDAgPyB2b2lkIDAgOiB0aWxlQm91bmRzLmJvdHRvbVJpZ2h0KSAhPT0gbnVsbCAmJiBfYiAhPT0gdm9pZCAwID8gX2IgOiAoMCwgdmVjXzEudmVjKShNYXRoLm1heCguLi4oX2QgPSAoX2MgPSBsYXllci5kYXRhKSA9PT0gbnVsbCB8fCBfYyA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2MubWFwKHJvdyA9PiByb3cubGVuZ3RoKSkgIT09IG51bGwgJiYgX2QgIT09IHZvaWQgMCA/IF9kIDogWzBdKSwgKF9mID0gKF9lID0gbGF5ZXIuZGF0YSkgPT09IG51bGwgfHwgX2UgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9lLmxlbmd0aCkgIT09IG51bGwgJiYgX2YgIT09IHZvaWQgMCA/IF9mIDogMCk7XG4gICAgICAgIGlmIChib3R0b21SaWdodC54IDw9IHRvcExlZnQueCB8fCBib3R0b21SaWdodC55IDw9IHRvcExlZnQueSkge1xuICAgICAgICAgICAgcmV0dXJuIFtdO1xuICAgICAgICB9XG4gICAgICAgIGNvbnN0IGJpdG1hcCA9IFtdO1xuICAgICAgICBmb3IgKGxldCB5ID0gdG9wTGVmdC55OyB5IDwgYm90dG9tUmlnaHQueTsgeSsrKSB7XG4gICAgICAgICAgICBjb25zdCByb3cgPSBbXTtcbiAgICAgICAgICAgIGZvciAobGV0IHggPSB0b3BMZWZ0Lng7IHggPCBib3R0b21SaWdodC54OyB4KyspIHtcbiAgICAgICAgICAgICAgICBjb25zdCB0aWxlRGF0YSA9IChfaCA9IChfZyA9IGxheWVyLmRhdGEpID09PSBudWxsIHx8IF9nID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfZ1t5XSkgPT09IG51bGwgfHwgX2ggPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9oW3hdO1xuICAgICAgICAgICAgICAgIGlmICh0aWxlRGF0YSA9PT0gdW5kZWZpbmVkIHx8IHRpbGVEYXRhID09PSAtMSkge1xuICAgICAgICAgICAgICAgICAgICByb3cucHVzaChmYWxzZSk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBjb25zdCB0aWxlID0gKF9qID0gbGF5ZXIudGlsZXMpID09PSBudWxsIHx8IF9qID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfalt0aWxlRGF0YV07XG4gICAgICAgICAgICAgICAgaWYgKCF0aWxlKSB7XG4gICAgICAgICAgICAgICAgICAgIHJvdy5wdXNoKGZhbHNlKTtcbiAgICAgICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGlmIChmaWVsZE5hbWUgJiYgIXRpbGVbZmllbGROYW1lXSkge1xuICAgICAgICAgICAgICAgICAgICByb3cucHVzaChmYWxzZSk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICByb3cucHVzaCh0cnVlKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGJpdG1hcC5wdXNoKHJvdyk7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuICgwLCBiaXRtYXBfZGVjb21wb3NlXzEuYml0bWFwVG9SZWN0YW5nbGVzKShiaXRtYXApO1xuICAgIH1cbiAgICAvKipcbiAgICAgKiBHZXQgdGhlIHRpbGUgYXQgYSBnaXZlbiBwb3NpdGlvbiBhbmQgaW4gdGhlIHNwZWNpZmllZCBsYXllclxuICAgICAqXG4gICAgICogSWYgbm8gbGF5ZXIgaXMgc3BlY2lmaWVkLCByZXR1cm4gYSBkaWN0aW9uYXJ5IG9mIGxheWVyIG5hbWVzIHRvIHRpbGVcbiAgICAgKiBkZWZpbml0aW9ucyAoaS5lLiByZXR1cm4gYWxsIGxheWVycylcbiAgICAgKlxuICAgICAqIElmIG5vIHRpbGUgZXhpc3RzIGF0IHRoaXMgcG9zaXRpb24sIHJldHVybiBudWxsXG4gICAgICovXG4gICAgZ2V0VGlsZUF0UG9zaXRpb24ocG9zaXRpb24sIGxheWVyTmFtZSkge1xuICAgICAgICBpZiAobGF5ZXJOYW1lKSB7XG4gICAgICAgICAgICByZXR1cm4gdGhpcy5nZXRUaWxlQXRQb3NpdGlvbkluTGF5ZXIocG9zaXRpb24sIGxheWVyTmFtZSk7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3QgcmVzdWx0ID0ge307XG4gICAgICAgIGZvciAoY29uc3QgbGF5ZXIgb2YgdGhpcy5vcHRpb25zLmxheWVycykge1xuICAgICAgICAgICAgcmVzdWx0W2xheWVyLm5hbWVdID0gdGhpcy5nZXRUaWxlQXRQb3NpdGlvbkluTGF5ZXIocG9zaXRpb24sIGxheWVyLm5hbWUpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuICAgIGdldFRpbGVBdFBvc2l0aW9uSW5MYXllcihwb3NpdGlvbiwgbGF5ZXJOYW1lKSB7XG4gICAgICAgIHZhciBfYSwgX2IsIF9jO1xuICAgICAgICBjb25zdCB0aWxlUG9zaXRpb24gPSB2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwocG9zaXRpb24sIDEgLyB0aGlzLm9wdGlvbnMudGlsZVNpemUpLCBNYXRoLmZsb29yKTtcbiAgICAgICAgY29uc3QgbGF5ZXIgPSB0aGlzLm9wdGlvbnMubGF5ZXJzLmZpbmQoKGwpID0+IGwubmFtZSA9PT0gbGF5ZXJOYW1lKTtcbiAgICAgICAgaWYgKCFsYXllcikge1xuICAgICAgICAgICAgcmV0dXJuIG51bGw7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3QgdGlsZURhdGEgPSAoX2IgPSAoX2EgPSBsYXllci5kYXRhKSA9PT0gbnVsbCB8fCBfYSA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2FbdGlsZVBvc2l0aW9uLnldKSA9PT0gbnVsbCB8fCBfYiA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2JbdGlsZVBvc2l0aW9uLnhdO1xuICAgICAgICBpZiAodGlsZURhdGEgPT09IHVuZGVmaW5lZCB8fCB0aWxlRGF0YSA9PT0gLTEpIHtcbiAgICAgICAgICAgIHJldHVybiBudWxsO1xuICAgICAgICB9XG4gICAgICAgIGlmIChsYXllci50aWxlcykge1xuICAgICAgICAgICAgcmV0dXJuIChfYyA9IGxheWVyLnRpbGVzW3RpbGVEYXRhXSkgIT09IG51bGwgJiYgX2MgIT09IHZvaWQgMCA/IF9jIDogbnVsbDtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gbnVsbDtcbiAgICB9XG4gICAgaGFzaFZlY3Rvcih2KSB7XG4gICAgICAgIHJldHVybiB2ZWNfMS52ZWMuc3RyKHYpO1xuICAgIH1cbiAgICBkcmF3KGNvbnRleHQsIHNjcmVlbiwgcG9zaXRpb24sIHNjYWxlKSB7XG4gICAgICAgIHZhciBfYSwgX2IsIF9jLCBfZDtcbiAgICAgICAgY29uc3QgYWJzb2x1dGVDaHVua1NpemUgPSB0aGlzLm9wdGlvbnMudGlsZVNpemUgKiB0aGlzLm9wdGlvbnMuY2h1bmtTaXplO1xuICAgICAgICBjb25zdCBjaHVua0JvcmRlciA9ICgwLCB2ZWNfMS52ZWMpKHRoaXMub3B0aW9ucy5jaHVua0JvcmRlcik7XG4gICAgICAgIC8vIE1heWJlIGNsYW1wIHNjYWxlXG4gICAgICAgIGxldCBhY3R1YWxTY2FsZSA9IHNjYWxlO1xuICAgICAgICBpZiAodGhpcy5vcHRpb25zLm1pblNjYWxlICYmIGFjdHVhbFNjYWxlIDwgdGhpcy5vcHRpb25zLm1pblNjYWxlKSB7XG4gICAgICAgICAgICBhY3R1YWxTY2FsZSA9IHRoaXMub3B0aW9ucy5taW5TY2FsZTtcbiAgICAgICAgfVxuICAgICAgICBpZiAodGhpcy5vcHRpb25zLm1heFNjYWxlICYmIGFjdHVhbFNjYWxlID4gdGhpcy5vcHRpb25zLm1heFNjYWxlKSB7XG4gICAgICAgICAgICBhY3R1YWxTY2FsZSA9IHRoaXMub3B0aW9ucy5tYXhTY2FsZTtcbiAgICAgICAgfVxuICAgICAgICAvLyBNYXliZSBjbGFtcCBwb3NpdGlvbiB0byBib3VuZHNcbiAgICAgICAgbGV0IGFjdHVhbFBvc2l0aW9uID0gKDAsIHZlY18xLnZlYykocG9zaXRpb24pO1xuICAgICAgICBpZiAodGhpcy5vcHRpb25zLmJvdW5kcyAmJiB0aGlzLm9wdGlvbnMuY2xhbXBQb3NpdGlvblRvQm91bmRzKSB7XG4gICAgICAgICAgICBjb25zdCB0aWxlU2l6ZVNjYWxlZCA9IHRoaXMub3B0aW9ucy50aWxlU2l6ZSAvIGFjdHVhbFNjYWxlO1xuICAgICAgICAgICAgY29uc3QgaGFsZlNjcmVlblNjYWxlZCA9IHZlY18xLnZlYy5tYXAodmVjXzEudmVjLm11bChzY3JlZW4sIDEgLyAoYWN0dWFsU2NhbGUgKiAyKSksIE1hdGguY2VpbCk7XG4gICAgICAgICAgICBjb25zdCBtaW5Qb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKHRoaXMub3B0aW9ucy5ib3VuZHMudG9wTGVmdC54ICogdGlsZVNpemVTY2FsZWQgKyBoYWxmU2NyZWVuU2NhbGVkLngsIHRoaXMub3B0aW9ucy5ib3VuZHMudG9wTGVmdC55ICogdGlsZVNpemVTY2FsZWQgKyBoYWxmU2NyZWVuU2NhbGVkLnkpO1xuICAgICAgICAgICAgY29uc3QgbWF4UG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh0aGlzLm9wdGlvbnMuYm91bmRzLmJvdHRvbVJpZ2h0LnggKiB0aWxlU2l6ZVNjYWxlZCAtIGhhbGZTY3JlZW5TY2FsZWQueCwgdGhpcy5vcHRpb25zLmJvdW5kcy5ib3R0b21SaWdodC55ICogdGlsZVNpemVTY2FsZWQgLSBoYWxmU2NyZWVuU2NhbGVkLnkpO1xuICAgICAgICAgICAgYWN0dWFsUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSgoMCwgdXRpbHNfMS5jbGFtcCkoYWN0dWFsUG9zaXRpb24ueCwgbWluUG9zaXRpb24ueCwgbWF4UG9zaXRpb24ueCksICgwLCB1dGlsc18xLmNsYW1wKShhY3R1YWxQb3NpdGlvbi55LCBtaW5Qb3NpdGlvbi55LCBtYXhQb3NpdGlvbi55KSk7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3Qgc2NyZWVuU2l6ZUluQ2h1bmtzID0gdmVjXzEudmVjLm1hcCh2ZWNfMS52ZWMubXVsKHNjcmVlbiwgMSAvIChhYnNvbHV0ZUNodW5rU2l6ZSAqIGFjdHVhbFNjYWxlKSksIE1hdGguY2VpbCk7XG4gICAgICAgIGNvbnN0IHNjcmVlbkNlbnRlckNodW5rID0gdmVjXzEudmVjLm1hcCh2ZWNfMS52ZWMubXVsKGFjdHVhbFBvc2l0aW9uLCAxIC8gYWJzb2x1dGVDaHVua1NpemUpLCBNYXRoLmZsb29yKTtcbiAgICAgICAgY29uc3QgdG9wTGVmdENodW5rID0gdmVjXzEudmVjLnN1Yih2ZWNfMS52ZWMuc3ViKHNjcmVlbkNlbnRlckNodW5rLCB2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwoc2NyZWVuU2l6ZUluQ2h1bmtzLCAwLjUpLCBNYXRoLmNlaWwpKSwgY2h1bmtCb3JkZXIpO1xuICAgICAgICBjb25zdCBib3R0b21SaWdodENodW5rID0gdmVjXzEudmVjLmFkZCh2ZWNfMS52ZWMuYWRkKHNjcmVlbkNlbnRlckNodW5rLCB2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwoc2NyZWVuU2l6ZUluQ2h1bmtzLCAwLjUpLCBNYXRoLmNlaWwpKSwgY2h1bmtCb3JkZXIpO1xuICAgICAgICBjb250ZXh0LnNhdmUoKTtcbiAgICAgICAgY29udGV4dC5zY2FsZShhY3R1YWxTY2FsZSwgYWN0dWFsU2NhbGUpO1xuICAgICAgICBjb250ZXh0LnRyYW5zbGF0ZSgtYWN0dWFsUG9zaXRpb24ueCArIHNjcmVlbi54IC8gKGFjdHVhbFNjYWxlICogMiksIC1hY3R1YWxQb3NpdGlvbi55ICsgc2NyZWVuLnkgLyAoYWN0dWFsU2NhbGUgKiAyKSk7XG4gICAgICAgIChfYiA9IChfYSA9IHRoaXMub3B0aW9ucykucHJlUmVuZGVyKSA9PT0gbnVsbCB8fCBfYiA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2IuY2FsbChfYSwgY29udGV4dCwgdGhpcywgc2NyZWVuLCBhY3R1YWxQb3NpdGlvbiwgYWN0dWFsU2NhbGUpO1xuICAgICAgICAvLyBSZW5kZXIgY2h1bmtzXG4gICAgICAgIGZvciAobGV0IHkgPSB0b3BMZWZ0Q2h1bmsueTsgeSA8IGJvdHRvbVJpZ2h0Q2h1bmsueTsgeSsrKSB7XG4gICAgICAgICAgICBmb3IgKGxldCB4ID0gdG9wTGVmdENodW5rLng7IHggPCBib3R0b21SaWdodENodW5rLng7IHgrKykge1xuICAgICAgICAgICAgICAgIGNvbnN0IGNodW5rUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh4LCB5KTtcbiAgICAgICAgICAgICAgICBjb25zdCBjaHVua0Fic29sdXRlUG9zaXRpb24gPSB2ZWNfMS52ZWMubXVsKGNodW5rUG9zaXRpb24sIGFic29sdXRlQ2h1bmtTaXplKTtcbiAgICAgICAgICAgICAgICAvLyBDaGVjayBpZiB3ZSBoYXZlIHRoaXMgY2h1bmsgaW4gdGhlIGNhY2hlXG4gICAgICAgICAgICAgICAgY29uc3QgY2h1bmtIYXNoID0gdGhpcy5oYXNoVmVjdG9yKGNodW5rUG9zaXRpb24pO1xuICAgICAgICAgICAgICAgIGlmICghdGhpcy5jaHVua0J1ZmZlci5oYXMoY2h1bmtIYXNoKSkge1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmNodW5rQnVmZmVyLnNldChjaHVua0hhc2gsIHRoaXMuZ2VuZXJhdGVDaHVuayhjaHVua1Bvc2l0aW9uLCBhYnNvbHV0ZUNodW5rU2l6ZSkpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBjb25zdCBjaHVuayA9IHRoaXMuY2h1bmtCdWZmZXIuZ2V0KGNodW5rSGFzaCk7XG4gICAgICAgICAgICAgICAgaWYgKGNodW5rKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnRleHQuZHJhd0ltYWdlKGNodW5rLmltYWdlLCBjaHVua0Fic29sdXRlUG9zaXRpb24ueCwgY2h1bmtBYnNvbHV0ZVBvc2l0aW9uLnkpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICAoX2QgPSAoX2MgPSB0aGlzLm9wdGlvbnMpLnBvc3RSZW5kZXIpID09PSBudWxsIHx8IF9kID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfZC5jYWxsKF9jLCBjb250ZXh0LCB0aGlzLCBzY3JlZW4sIGFjdHVhbFBvc2l0aW9uLCBhY3R1YWxTY2FsZSk7XG4gICAgICAgIC8vIFJlbmRlciBkZWJ1ZyBoZWxwZXJzXG4gICAgICAgIGlmICh0aGlzLm9wdGlvbnMuZGVidWcuc2hvd1RpbGVCb3JkZXJzKSB7XG4gICAgICAgICAgICBjb25zdCB0b3BMZWZ0VGlsZSA9IHZlY18xLnZlYy5tdWwodmVjXzEudmVjLnN1YihzY3JlZW5DZW50ZXJDaHVuaywgdmVjXzEudmVjLmFkZCh2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwoc2NyZWVuU2l6ZUluQ2h1bmtzLCAwLjUpLCBNYXRoLmNlaWwpLCAoMCwgdmVjXzEudmVjKSgxKSkpLCB0aGlzLm9wdGlvbnMuY2h1bmtTaXplKTtcbiAgICAgICAgICAgIGNvbnN0IGJvdHRvbVJpZ2h0VGlsZSA9IHZlY18xLnZlYy5tdWwodmVjXzEudmVjLmFkZChzY3JlZW5DZW50ZXJDaHVuaywgdmVjXzEudmVjLmFkZCh2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwoc2NyZWVuU2l6ZUluQ2h1bmtzLCAwLjUpLCBNYXRoLmNlaWwpLCAoMCwgdmVjXzEudmVjKSgxKSkpLCB0aGlzLm9wdGlvbnMuY2h1bmtTaXplKTtcbiAgICAgICAgICAgIGZvciAobGV0IHkgPSB0b3BMZWZ0VGlsZS55OyB5IDwgYm90dG9tUmlnaHRUaWxlLnk7IHkrKykge1xuICAgICAgICAgICAgICAgIHRoaXMuZHJhd0xpbmUoY29udGV4dCwgKDAsIHZlY18xLnZlYykoYWN0dWFsUG9zaXRpb24ueCAtIHNjcmVlbi54IC8gKGFjdHVhbFNjYWxlICogMiksIHkgKiB0aGlzLm9wdGlvbnMudGlsZVNpemUpLCAoMCwgdmVjXzEudmVjKShhY3R1YWxQb3NpdGlvbi54ICsgc2NyZWVuLnggLyAoYWN0dWFsU2NhbGUgKiAyKSwgeSAqIHRoaXMub3B0aW9ucy50aWxlU2l6ZSksIFRpbGVNYXAuREVCVUdfVElMRV9CT1JERVJfQ09MT1VSLCBUaWxlTWFwLkRFQlVHX1RJTEVfQk9SREVSX0xJTkVfV0lEVEgpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZm9yIChsZXQgeCA9IHRvcExlZnRUaWxlLng7IHggPCBib3R0b21SaWdodFRpbGUueDsgeCsrKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5kcmF3TGluZShjb250ZXh0LCAoMCwgdmVjXzEudmVjKSh4ICogdGhpcy5vcHRpb25zLnRpbGVTaXplLCBhY3R1YWxQb3NpdGlvbi55IC0gc2NyZWVuLnkgLyAoYWN0dWFsU2NhbGUgKiAyKSksICgwLCB2ZWNfMS52ZWMpKHggKiB0aGlzLm9wdGlvbnMudGlsZVNpemUsIGFjdHVhbFBvc2l0aW9uLnkgKyBzY3JlZW4ueSAvIChhY3R1YWxTY2FsZSAqIDIpKSwgVGlsZU1hcC5ERUJVR19USUxFX0JPUkRFUl9DT0xPVVIsIFRpbGVNYXAuREVCVUdfVElMRV9CT1JERVJfTElORV9XSURUSCk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5kZWJ1Zy5zaG93Q2h1bmtCb3JkZXJzKSB7XG4gICAgICAgICAgICBmb3IgKGxldCB5ID0gdG9wTGVmdENodW5rLnk7IHkgPCBib3R0b21SaWdodENodW5rLnk7IHkrKykge1xuICAgICAgICAgICAgICAgIHRoaXMuZHJhd0xpbmUoY29udGV4dCwgKDAsIHZlY18xLnZlYykoYWN0dWFsUG9zaXRpb24ueCAtIHNjcmVlbi54IC8gKGFjdHVhbFNjYWxlICogMiksIHkgKiBhYnNvbHV0ZUNodW5rU2l6ZSksICgwLCB2ZWNfMS52ZWMpKGFjdHVhbFBvc2l0aW9uLnggKyBzY3JlZW4ueCAvIChhY3R1YWxTY2FsZSAqIDIpLCB5ICogYWJzb2x1dGVDaHVua1NpemUpLCBUaWxlTWFwLkRFQlVHX0NIVU5LX0JPUkRFUl9DT0xPVVIsIFRpbGVNYXAuREVCVUdfQ0hVTktfQk9SREVSX0xJTkVfV0lEVEgpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZm9yIChsZXQgeCA9IHRvcExlZnRDaHVuay54OyB4IDwgYm90dG9tUmlnaHRDaHVuay54OyB4KyspIHtcbiAgICAgICAgICAgICAgICB0aGlzLmRyYXdMaW5lKGNvbnRleHQsICgwLCB2ZWNfMS52ZWMpKHggKiBhYnNvbHV0ZUNodW5rU2l6ZSwgYWN0dWFsUG9zaXRpb24ueSAtIHNjcmVlbi55IC8gKGFjdHVhbFNjYWxlICogMikpLCAoMCwgdmVjXzEudmVjKSh4ICogYWJzb2x1dGVDaHVua1NpemUsIGFjdHVhbFBvc2l0aW9uLnkgKyBzY3JlZW4ueSAvIChhY3R1YWxTY2FsZSAqIDIpKSwgVGlsZU1hcC5ERUJVR19DSFVOS19CT1JERVJfQ09MT1VSLCBUaWxlTWFwLkRFQlVHX0NIVU5LX0JPUkRFUl9MSU5FX1dJRFRIKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBpZiAodGhpcy5vcHRpb25zLmRlYnVnLnNob3dDaHVua0xhYmVscykge1xuICAgICAgICAgICAgY29udGV4dC5zYXZlKCk7XG4gICAgICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFRpbGVNYXAuREVCVUdfQ0hVTktfTEFCRUxfQ09MT1VSO1xuICAgICAgICAgICAgY29udGV4dC5mb250ID0gVGlsZU1hcC5ERUJVR19DSFVOS19MQUJFTF9GT05UO1xuICAgICAgICAgICAgY29udGV4dC50ZXh0QmFzZWxpbmUgPSAnbWlkZGxlJztcbiAgICAgICAgICAgIGNvbnRleHQudGV4dEFsaWduID0gJ2NlbnRlcic7XG4gICAgICAgICAgICBmb3IgKGxldCB5ID0gdG9wTGVmdENodW5rLnk7IHkgPCBib3R0b21SaWdodENodW5rLnk7IHkrKykge1xuICAgICAgICAgICAgICAgIGZvciAobGV0IHggPSB0b3BMZWZ0Q2h1bmsueDsgeCA8IGJvdHRvbVJpZ2h0Q2h1bmsueDsgeCsrKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnRleHQuZmlsbFRleHQoYCR7eH0sICR7eX1gLCB4ICogYWJzb2x1dGVDaHVua1NpemUgKyBhYnNvbHV0ZUNodW5rU2l6ZSAvIDIsIHkgKiBhYnNvbHV0ZUNodW5rU2l6ZSArIGFic29sdXRlQ2h1bmtTaXplIC8gMik7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgY29udGV4dC5yZXN0b3JlKCk7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5kZWJ1Zy5zaG93T3JpZ2luICYmXG4gICAgICAgICAgICBwb2ludEluUmVjdGFuZ2xlKCgwLCB2ZWNfMS52ZWMpKDAsIDApLCB0b3BMZWZ0Q2h1bmssIGJvdHRvbVJpZ2h0Q2h1bmspKSB7XG4gICAgICAgICAgICB0aGlzLmRyYXdDcm9zcyhjb250ZXh0LCAoMCwgdmVjXzEudmVjKSgwLCAwKSwgVGlsZU1hcC5ERUJVR19PUklHSU5fQ09MT1VSLCBUaWxlTWFwLkRFQlVHX09SSUdJTl9MSU5FX1dJRFRILCBUaWxlTWFwLkRFQlVHX09SSUdJTl9TSVpFKTtcbiAgICAgICAgfVxuICAgICAgICBjb250ZXh0LnJlc3RvcmUoKTtcbiAgICB9XG4gICAgZ2VuZXJhdGVDaHVuayhjaHVua1Bvc2l0aW9uLCBhYnNvbHV0ZUNodW5rU2l6ZSkge1xuICAgICAgICB2YXIgX2EsIF9iLCBfYywgX2QsIF9lLCBfZiwgX2csIF9oLCBfaiwgX2ssIF9sLCBfbTtcbiAgICAgICAgY29uc3QgY2h1bmtDYW52YXMgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KCdjYW52YXMnKTtcbiAgICAgICAgY29uc3QgY2h1bmtDb250ZXh0ID0gY2h1bmtDYW52YXMuZ2V0Q29udGV4dCgnMmQnKTtcbiAgICAgICAgY2h1bmtDYW52YXMud2lkdGggPSBhYnNvbHV0ZUNodW5rU2l6ZTtcbiAgICAgICAgY2h1bmtDYW52YXMuaGVpZ2h0ID0gYWJzb2x1dGVDaHVua1NpemU7XG4gICAgICAgIGxldCBjaHVuayA9IHtcbiAgICAgICAgICAgIGNodW5rUG9zaXRpb24sXG4gICAgICAgICAgICBpbWFnZTogY2h1bmtDYW52YXMsXG4gICAgICAgIH07XG4gICAgICAgIGNvbnN0IHRvcExlZnRUaWxlID0gdmVjXzEudmVjLm11bChjaHVua1Bvc2l0aW9uLCB0aGlzLm9wdGlvbnMuY2h1bmtTaXplKTtcbiAgICAgICAgY29uc3QgYm90dG9tUmlnaHRUaWxlID0gdmVjXzEudmVjLmFkZCh0b3BMZWZ0VGlsZSwgKDAsIHZlY18xLnZlYykodGhpcy5vcHRpb25zLmNodW5rU2l6ZSAtIDEpKTtcbiAgICAgICAgY29uc3QgYm91bmRzVG9wTGVmdCA9IChfYiA9IChfYSA9IHRoaXMub3B0aW9ucy5ib3VuZHMpID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS50b3BMZWZ0KSAhPT0gbnVsbCAmJiBfYiAhPT0gdm9pZCAwID8gX2IgOiAoMCwgdmVjXzEudmVjKSgwKTtcbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5wcmVHZW5lcmF0ZUNodW5rKSB7XG4gICAgICAgICAgICBjb25zdCByZXN1bHQgPSB0aGlzLm9wdGlvbnMucHJlR2VuZXJhdGVDaHVuayhjaHVua0NvbnRleHQsIHRoaXMsIHtcbiAgICAgICAgICAgICAgICB0b3BMZWZ0OiB0b3BMZWZ0VGlsZSxcbiAgICAgICAgICAgICAgICBib3R0b21SaWdodDogYm90dG9tUmlnaHRUaWxlLFxuICAgICAgICAgICAgfSwgY2h1bmtQb3NpdGlvbik7XG4gICAgICAgICAgICBpZiAoQXJyYXkuaXNBcnJheShyZXN1bHQpKSB7XG4gICAgICAgICAgICAgICAgaWYgKCFyZXN1bHRbMV0pIHtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNodW5rO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICAvLyBEZWZhdWx0IGdlbmVyYXRpb24sIHJlbmRlciB0aWxlcyBmcm9tIHRpbGVtYXAgZGF0YVxuICAgICAgICBmb3IgKGNvbnN0IGxheWVyIG9mIHRoaXMub3B0aW9ucy5sYXllcnMpIHtcbiAgICAgICAgICAgIGNodW5rQ29udGV4dC5zYXZlKCk7XG4gICAgICAgICAgICBjaHVua0NvbnRleHQuZ2xvYmFsQWxwaGEgPSAoX2MgPSBsYXllci5vcGFjaXR5KSAhPT0gbnVsbCAmJiBfYyAhPT0gdm9pZCAwID8gX2MgOiAxO1xuICAgICAgICAgICAgY29uc3QgYWxpZ25tZW50ID0gKF9kID0gbGF5ZXIuYWxpZ25tZW50KSAhPT0gbnVsbCAmJiBfZCAhPT0gdm9pZCAwID8gX2QgOiBUaWxlQWxpZ25tZW50LkNlbnRlcjtcbiAgICAgICAgICAgIGZvciAobGV0IHkgPSB0b3BMZWZ0VGlsZS55OyB5IDw9IGJvdHRvbVJpZ2h0VGlsZS55OyB5KyspIHtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCB4ID0gdG9wTGVmdFRpbGUueDsgeCA8PSBib3R0b21SaWdodFRpbGUueDsgeCsrKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IHRpbGVQb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKHgsIHkpO1xuICAgICAgICAgICAgICAgICAgICAoX2UgPSBsYXllci5wcmVSZW5kZXJUaWxlKSA9PT0gbnVsbCB8fCBfZSA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2UuY2FsbChsYXllciwgY2h1bmtDb250ZXh0LCB0aGlzLCBsYXllciwgY2h1bmtQb3NpdGlvbiwgdGlsZVBvc2l0aW9uKTtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgdGlsZURhdGFQb3NpdGlvbiA9IHZlY18xLnZlYy5zdWIodGlsZVBvc2l0aW9uLCBib3VuZHNUb3BMZWZ0KTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKHRpbGVEYXRhUG9zaXRpb24ueCA8IDAgfHwgdGlsZURhdGFQb3NpdGlvbi55IDwgMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgY29uc3QgdGlsZURhdGEgPSAoX2cgPSAoX2YgPSBsYXllci5kYXRhKSA9PT0gbnVsbCB8fCBfZiA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2ZbdGlsZURhdGFQb3NpdGlvbi55XSkgPT09IG51bGwgfHwgX2cgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9nW3RpbGVEYXRhUG9zaXRpb24ueF07XG4gICAgICAgICAgICAgICAgICAgIGlmICh0aWxlRGF0YSA9PT0gdW5kZWZpbmVkIHx8IHRpbGVEYXRhID09PSAtMSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgY29uc3QgdGlsZUltYWdlID0gKF9qID0gKF9oID0gbGF5ZXIudGlsZXMpID09PSBudWxsIHx8IF9oID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfaFt0aWxlRGF0YV0pID09PSBudWxsIHx8IF9qID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfai5pbWFnZTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKCF0aWxlSW1hZ2UpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IHRpbGVBYnNvbHV0ZVBvc2l0aW9uID0gdmVjXzEudmVjLnN1Yih2ZWNfMS52ZWMubXVsKHRpbGVQb3NpdGlvbiwgdGhpcy5vcHRpb25zLnRpbGVTaXplKSwgdmVjXzEudmVjLm11bChjaHVua1Bvc2l0aW9uLCBhYnNvbHV0ZUNodW5rU2l6ZSkpO1xuICAgICAgICAgICAgICAgICAgICAvLyBUaWxlIGNsaXBwaW5nXG4gICAgICAgICAgICAgICAgICAgIGlmIChsYXllci5jbGlwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBjaHVua0NvbnRleHQuc2F2ZSgpO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2h1bmtDb250ZXh0LmJlZ2luUGF0aCgpO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2h1bmtDb250ZXh0LnJlY3QodGlsZUFic29sdXRlUG9zaXRpb24ueCwgdGlsZUFic29sdXRlUG9zaXRpb24ueSwgdGhpcy5vcHRpb25zLnRpbGVTaXplLCB0aGlzLm9wdGlvbnMudGlsZVNpemUpO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2h1bmtDb250ZXh0LmNsaXAoKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAvLyBUaWxlIGFsaWdubWVudFxuICAgICAgICAgICAgICAgICAgICBsZXQgdGlsZUltYWdlQWJzb2x1dGVQb3NpdGlvbjtcbiAgICAgICAgICAgICAgICAgICAgc3dpdGNoIChhbGlnbm1lbnQpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNhc2UgVGlsZUFsaWdubWVudC5Ub3BMZWZ0OlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh0aWxlQWJzb2x1dGVQb3NpdGlvbik7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFRpbGVBbGlnbm1lbnQuVG9wOlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSgodGlsZUFic29sdXRlUG9zaXRpb24ueCArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAvIDIpIC0gdGlsZUltYWdlLndpZHRoIC8gMiwgdGlsZUFic29sdXRlUG9zaXRpb24ueSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFRpbGVBbGlnbm1lbnQuVG9wUmlnaHQ6XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdGlsZUltYWdlQWJzb2x1dGVQb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKHRpbGVBYnNvbHV0ZVBvc2l0aW9uLnggKyB0aGlzLm9wdGlvbnMudGlsZVNpemUgLSB0aWxlSW1hZ2Uud2lkdGgsIHRpbGVBYnNvbHV0ZVBvc2l0aW9uLnkpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2FzZSBUaWxlQWxpZ25tZW50LkxlZnQ6XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdGlsZUltYWdlQWJzb2x1dGVQb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKHRpbGVBYnNvbHV0ZVBvc2l0aW9uLngsICh0aWxlQWJzb2x1dGVQb3NpdGlvbi55ICsgdGhpcy5vcHRpb25zLnRpbGVTaXplIC8gMikgLSB0aWxlSW1hZ2UuaGVpZ2h0IC8gMik7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFRpbGVBbGlnbm1lbnQuQ2VudGVyOlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSgodGlsZUFic29sdXRlUG9zaXRpb24ueCArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAvIDIpIC0gdGlsZUltYWdlLndpZHRoIC8gMiwgKHRpbGVBYnNvbHV0ZVBvc2l0aW9uLnkgKyB0aGlzLm9wdGlvbnMudGlsZVNpemUgLyAyKSAtIHRpbGVJbWFnZS5oZWlnaHQgLyAyKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICAgICAgICAgIGNhc2UgVGlsZUFsaWdubWVudC5SaWdodDpcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB0aWxlSW1hZ2VBYnNvbHV0ZVBvc2l0aW9uID0gKDAsIHZlY18xLnZlYykodGlsZUFic29sdXRlUG9zaXRpb24ueCArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAtIHRpbGVJbWFnZS53aWR0aCwgKHRpbGVBYnNvbHV0ZVBvc2l0aW9uLnkgKyB0aGlzLm9wdGlvbnMudGlsZVNpemUgLyAyKSAtIHRpbGVJbWFnZS5oZWlnaHQgLyAyKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICAgICAgICAgIGNhc2UgVGlsZUFsaWdubWVudC5Cb3R0b21MZWZ0OlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh0aWxlQWJzb2x1dGVQb3NpdGlvbi54LCB0aWxlQWJzb2x1dGVQb3NpdGlvbi55ICsgdGhpcy5vcHRpb25zLnRpbGVTaXplIC0gdGlsZUltYWdlLmhlaWdodCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFRpbGVBbGlnbm1lbnQuQm90dG9tOlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSgodGlsZUFic29sdXRlUG9zaXRpb24ueCArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAvIDIpIC0gdGlsZUltYWdlLndpZHRoIC8gMiwgdGlsZUFic29sdXRlUG9zaXRpb24ueSArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAtIHRpbGVJbWFnZS5oZWlnaHQpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2FzZSBUaWxlQWxpZ25tZW50LkJvdHRvbVJpZ2h0OlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh0aWxlQWJzb2x1dGVQb3NpdGlvbi54ICsgdGhpcy5vcHRpb25zLnRpbGVTaXplIC0gdGlsZUltYWdlLndpZHRoLCB0aWxlQWJzb2x1dGVQb3NpdGlvbi55ICsgdGhpcy5vcHRpb25zLnRpbGVTaXplIC0gdGlsZUltYWdlLmhlaWdodCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgY2h1bmtDb250ZXh0LmRyYXdJbWFnZSh0aWxlSW1hZ2UsIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24ueCwgdGlsZUltYWdlQWJzb2x1dGVQb3NpdGlvbi55KTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGxheWVyLmNsaXApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNodW5rQ29udGV4dC5yZXN0b3JlKCk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgKF9rID0gbGF5ZXIucG9zdFJlbmRlclRpbGUpID09PSBudWxsIHx8IF9rID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfay5jYWxsKGxheWVyLCBjaHVua0NhbnZhcywgY2h1bmtDb250ZXh0LCB0aGlzLCBsYXllciwgY2h1bmtQb3NpdGlvbiwgdGlsZVBvc2l0aW9uKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBjaHVua0NvbnRleHQucmVzdG9yZSgpO1xuICAgICAgICB9XG4gICAgICAgIChfbSA9IChfbCA9IHRoaXMub3B0aW9ucykucG9zdEdlbmVyYXRlQ2h1bmspID09PSBudWxsIHx8IF9tID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfbS5jYWxsKF9sLCBjaHVua0NhbnZhcywgY2h1bmtDb250ZXh0LCB0aGlzLCB7XG4gICAgICAgICAgICB0b3BMZWZ0OiB0b3BMZWZ0VGlsZSxcbiAgICAgICAgICAgIGJvdHRvbVJpZ2h0OiBib3R0b21SaWdodFRpbGUsXG4gICAgICAgIH0sIGNodW5rUG9zaXRpb24pO1xuICAgICAgICByZXR1cm4gY2h1bms7XG4gICAgfVxuICAgIGRyYXdMaW5lKGNvbnRleHQsIHN0YXJ0LCBlbmQsIGNvbG91ciwgbGluZVdpZHRoKSB7XG4gICAgICAgIGNvbnRleHQuc2F2ZSgpO1xuICAgICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IGxpbmVXaWR0aDtcbiAgICAgICAgY29udGV4dC5zdHJva2VTdHlsZSA9IGNvbG91cjtcbiAgICAgICAgY29udGV4dC5iZWdpblBhdGgoKTtcbiAgICAgICAgY29udGV4dC5tb3ZlVG8oc3RhcnQueCwgc3RhcnQueSk7XG4gICAgICAgIGNvbnRleHQubGluZVRvKGVuZC54LCBlbmQueSk7XG4gICAgICAgIGNvbnRleHQuc3Ryb2tlKCk7XG4gICAgICAgIGNvbnRleHQucmVzdG9yZSgpO1xuICAgIH1cbiAgICBkcmF3Q3Jvc3MoY29udGV4dCwgcG9zaXRpb24sIGNvbG91ciwgbGluZVdpZHRoLCBzaXplKSB7XG4gICAgICAgIGNvbnRleHQuc2F2ZSgpO1xuICAgICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IGxpbmVXaWR0aDtcbiAgICAgICAgY29uc3QgaGFsZlNpemUgPSBNYXRoLmNlaWwoc2l6ZSAvIDIpO1xuICAgICAgICBjb250ZXh0LnN0cm9rZVN0eWxlID0gY29sb3VyO1xuICAgICAgICBjb250ZXh0LmJlZ2luUGF0aCgpO1xuICAgICAgICBjb250ZXh0Lm1vdmVUbyhwb3NpdGlvbi54IC0gaGFsZlNpemUsIHBvc2l0aW9uLnkgLSBoYWxmU2l6ZSk7XG4gICAgICAgIGNvbnRleHQubGluZVRvKHBvc2l0aW9uLnggKyBoYWxmU2l6ZSwgcG9zaXRpb24ueSArIGhhbGZTaXplKTtcbiAgICAgICAgY29udGV4dC5tb3ZlVG8ocG9zaXRpb24ueCAtIGhhbGZTaXplLCBwb3NpdGlvbi55ICsgaGFsZlNpemUpO1xuICAgICAgICBjb250ZXh0LmxpbmVUbyhwb3NpdGlvbi54ICsgaGFsZlNpemUsIHBvc2l0aW9uLnkgLSBoYWxmU2l6ZSk7XG4gICAgICAgIGNvbnRleHQuc3Ryb2tlKCk7XG4gICAgICAgIGNvbnRleHQucmVzdG9yZSgpO1xuICAgIH1cbn1cbmV4cG9ydHMuVGlsZU1hcCA9IFRpbGVNYXA7XG5UaWxlTWFwLkRFRkFVTFRfT1BUSU9OUyA9IHtcbiAgICBjbGFtcFBvc2l0aW9uVG9Cb3VuZHM6IHRydWUsXG4gICAgdGlsZVNpemU6IDE2LFxuICAgIGxheWVyczogW1xuICAgICAgICB7XG4gICAgICAgICAgICBuYW1lOiAnZGVmYXVsdCcsXG4gICAgICAgIH0sXG4gICAgXSxcbiAgICBjaHVua1NpemU6IDgsXG4gICAgY2h1bmtCb3JkZXI6IDEsXG4gICAgY2h1bmtCdWZmZXJNYXhTaXplOiA2NCxcbn07XG5UaWxlTWFwLkRFQlVHX09SSUdJTl9DT0xPVVIgPSAnY3lhbic7XG5UaWxlTWFwLkRFQlVHX09SSUdJTl9MSU5FX1dJRFRIID0gMjtcblRpbGVNYXAuREVCVUdfT1JJR0lOX1NJWkUgPSAxMDtcblRpbGVNYXAuREVCVUdfQ0hVTktfQk9SREVSX0NPTE9VUiA9ICd5ZWxsb3cnO1xuVGlsZU1hcC5ERUJVR19DSFVOS19CT1JERVJfTElORV9XSURUSCA9IDI7XG5UaWxlTWFwLkRFQlVHX0NIVU5LX0xBQkVMX0NPTE9VUiA9ICd3aGl0ZSc7XG5UaWxlTWFwLkRFQlVHX0NIVU5LX0xBQkVMX0ZPTlQgPSAnMTJweCBtb25vc3BhY2UnO1xuVGlsZU1hcC5ERUJVR19USUxFX0JPUkRFUl9DT0xPVVIgPSAnb3JhbmdlJztcblRpbGVNYXAuREVCVUdfVElMRV9CT1JERVJfTElORV9XSURUSCA9IDE7XG4vKipcbiAqIENvbnRlbnQgTWFuYWdlciBQcm9jZXNzb3Igd3JhcHBlciB3aGljaCBjb252ZXJ0cyBUaWxlTWFwT3B0aW9uc0RhdGEgaW50b1xuICogVGlsZU1hcE9wdGlvbnNcbiAqXG4gKiBAc2VlIGh0dHBzOi8vd3d3Lm5wbWpzLmNvbS9wYWNrYWdlL0BiYXNlbWVudHVuaXZlcnNlL2NvbnRlbnQtbWFuYWdlclxuICovXG5hc3luYyBmdW5jdGlvbiB0aWxlTWFwT3B0aW9uc0NvbnRlbnRQcm9jZXNzb3IoY29udGVudCwgZGF0YSwgb3B0aW9ucykge1xuICAgIGNvbnN0IGdldEltYWdlRnJvbUNvbnRlbnQgPSAobmFtZSkgPT4ge1xuICAgICAgICB2YXIgX2E7XG4gICAgICAgIGNvbnN0IGltYWdlID0gKF9hID0gY29udGVudFtuYW1lXSkgPT09IG51bGwgfHwgX2EgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9hLmNvbnRlbnQ7XG4gICAgICAgIGlmICghaW1hZ2UpIHtcbiAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihgSW1hZ2UgJyR7bmFtZX0nIG5vdCBmb3VuZGApO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBpbWFnZTtcbiAgICB9O1xuICAgIGNvbnN0IHJlc3VsdCA9IGRhdGE7XG4gICAgaWYgKHJlc3VsdC5sYXllcnMpIHtcbiAgICAgICAgZm9yIChjb25zdCBbaSwgbGF5ZXJdIG9mIHJlc3VsdC5sYXllcnMuZW50cmllcygpKSB7XG4gICAgICAgICAgICAvLyBSZXBsYWNlIGltYWdlTmFtZSB3aXRoIGltYWdlIGluIHRoZSB0aWxlIGRlZmluaXRpb25zIGFycmF5XG4gICAgICAgICAgICBpZiAobGF5ZXIudGlsZXMpIHtcbiAgICAgICAgICAgICAgICBmb3IgKGNvbnN0IFtqLCB0aWxlXSBvZiBsYXllci50aWxlcy5lbnRyaWVzKCkpIHtcbiAgICAgICAgICAgICAgICAgICAgcmVzdWx0LmxheWVyc1tpXS50aWxlc1tqXS5pbWFnZSA9IGdldEltYWdlRnJvbUNvbnRlbnQodGlsZS5pbWFnZU5hbWUpO1xuICAgICAgICAgICAgICAgICAgICBkZWxldGUgcmVzdWx0LmxheWVyc1tpXS50aWxlc1tqXS5pbWFnZU5hbWU7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgLy8gRGVjb21wcmVzcyBsYXllciBkYXRhXG4gICAgICAgICAgICBpZiAoKG9wdGlvbnMgPT09IG51bGwgfHwgb3B0aW9ucyA9PT0gdm9pZCAwID8gdm9pZCAwIDogb3B0aW9ucy5kZWNvbXByZXNzRGF0YSkgJiYgbGF5ZXIuZGF0YSAmJiBsYXllci53aWR0aCkge1xuICAgICAgICAgICAgICAgIHJlc3VsdC5sYXllcnNbaV0uZGF0YSA9ICgwLCB1dGlsc18xLmNodW5rKSgoMCwgZGVjb2RlXzEuZGVmYXVsdCkobGF5ZXIuZGF0YSksIGxheWVyLndpZHRoKTtcbiAgICAgICAgICAgICAgICBkZWxldGUgcmVzdWx0LmxheWVyc1tpXS53aWR0aDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cbiAgICAvLyBAdHMtaWdub3JlXG4gICAgZGF0YS5jb250ZW50ID0gcmVzdWx0O1xufVxuZXhwb3J0cy50aWxlTWFwT3B0aW9uc0NvbnRlbnRQcm9jZXNzb3IgPSB0aWxlTWFwT3B0aW9uc0NvbnRlbnRQcm9jZXNzb3I7XG4vLyMgc291cmNlTWFwcGluZ1VSTD1kYXRhOmFwcGxpY2F0aW9uL2pzb247YmFzZTY0LGV5SjJaWEp6YVc5dUlqb3pMQ0ptYVd4bElqb2lhVzVrWlhndWFuTWlMQ0p6YjNWeVkyVlNiMjkwSWpvaUlpd2ljMjkxY21ObGN5STZXeUl1TGk5cGJtUmxlQzUwY3lKZExDSnVZVzFsY3lJNlcxMHNJbTFoY0hCcGJtZHpJam9pT3pzN096czdRVUZCUVN4eFEwRkJhVU03UVVGRGFrTXNOa1JCUVhGRE8wRkJRM0pETEN0RFFVRTBRenRCUVVNMVF5eHRSRUZCZFVRN1FVRkRka1FzZVVSQlFXMUZPMEZCTkZGdVJTeEpRVUZaTEdGQlZWZzdRVUZXUkN4WFFVRlpMR0ZCUVdFN1NVRkRka0lzZFVSQlFWY3NRMEZCUVR0SlFVTllMQ3REUVVGSExFTkJRVUU3U1VGRFNDeDVSRUZCVVN4RFFVRkJPMGxCUTFJc2FVUkJRVWtzUTBGQlFUdEpRVU5LTEhGRVFVRk5MRU5CUVVFN1NVRkRUaXh0UkVGQlN5eERRVUZCTzBsQlEwd3NOa1JCUVZVc1EwRkJRVHRKUVVOV0xIRkVRVUZOTEVOQlFVRTdTVUZEVGl3clJFRkJWeXhEUVVGQk8wRkJRMklzUTBGQlF5eEZRVlpYTEdGQlFXRXNSMEZCWWl4eFFrRkJZU3hMUVVGaUxIRkNRVUZoTEZGQlZYaENPMEZCWlVRc1UwRkJVeXhuUWtGQlowSXNRMEZEZGtJc1MwRkJWU3hGUVVOV0xFOUJRVmtzUlVGRFdpeFhRVUZuUWp0SlFVVm9RaXhQUVVGUExFTkJRMHdzUzBGQlN5eERRVUZETEVOQlFVTXNTVUZCU1N4UFFVRlBMRU5CUVVNc1EwRkJRenRSUVVOd1FpeExRVUZMTEVOQlFVTXNRMEZCUXl4SlFVRkpMRTlCUVU4c1EwRkJReXhEUVVGRE8xRkJRM0JDTEV0QlFVc3NRMEZCUXl4RFFVRkRMRWRCUVVjc1YwRkJWeXhEUVVGRExFTkJRVU03VVVGRGRrSXNTMEZCU3l4RFFVRkRMRU5CUVVNc1IwRkJSeXhYUVVGWExFTkJRVU1zUTBGQlF5eERRVU40UWl4RFFVRkRPMEZCUTBvc1EwRkJRenRCUVVWRUxFMUJRV0VzVDBGQlR6dEpRV2xEYkVJc1dVRkJiVUlzVDBGQmIwTTdVVUZEY2tRc1RVRkJUU3hoUVVGaExFZEJRVWNzVFVGQlRTeERRVUZETEUxQlFVMHNRMEZEYWtNc1JVRkJSU3hGUVVOR0xFOUJRVThzUTBGQlF5eGxRVUZsTEVWQlEzWkNMRTlCUVU4c1lVRkJVQ3hQUVVGUExHTkJRVkFzVDBGQlR5eEhRVUZKTEVWQlFVVXNRMEZEWkN4RFFVRkRPMUZCUlVZc1NVRkJTU3hEUVVGRExHRkJRV0VzUTBGQlF5eExRVUZMTEVsQlFVa3NZVUZCWVN4RFFVRkRMRXRCUVVzc1MwRkJTeXhKUVVGSkxFVkJRVVU3V1VGRGVFUXNZVUZCWVN4RFFVRkRMRXRCUVVzc1IwRkJSenRuUWtGRGNFSXNWVUZCVlN4RlFVRkZMRU5CUVVNc1EwRkJReXhoUVVGaExFTkJRVU1zUzBGQlN6dG5Ra0ZEYWtNc1owSkJRV2RDTEVWQlFVVXNRMEZCUXl4RFFVRkRMR0ZCUVdFc1EwRkJReXhMUVVGTE8yZENRVU4yUXl4bFFVRmxMRVZCUVVVc1EwRkJReXhEUVVGRExHRkJRV0VzUTBGQlF5eExRVUZMTzJkQ1FVTjBReXhsUVVGbExFVkJRVVVzUTBGQlF5eERRVUZETEdGQlFXRXNRMEZCUXl4TFFVRkxPMkZCUTNaRExFTkJRVU03VTBGRFNEdFJRVVZFTEVsQlFVa3NRMEZCUXl4UFFVRlBMRWRCUVVjc1lVRkJiME1zUTBGQlF6dFJRVVZ3UkN4SlFVRkpMRU5CUVVNc1YwRkJWeXhIUVVGSExFbEJRVWtzWjBKQlFVMHNRMEZCUXl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExHdENRVUZyUWl4RFFVRkRMRU5CUVVNN1NVRkRha1VzUTBGQlF6dEpRVVZFT3pzN096czdPenM3VDBGVFJ6dEpRVU5KTEd0Q1FVRnJRaXhEUVVOMlFpeFRRVUZwUWl4RlFVTnFRaXhUUVVGdFF5eEZRVU51UXl4VlFVRnRRanM3VVVGRmJrSXNUVUZCVFN4TFFVRkxMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eE5RVUZOTEVOQlFVTXNTVUZCU1N4RFFVRkRMRU5CUVVNc1EwRkJReXhGUVVGRkxFVkJRVVVzUTBGQlF5eERRVUZETEVOQlFVTXNTVUZCU1N4TFFVRkxMRk5CUVZNc1EwRkJReXhEUVVGRE8xRkJRM0JGTEVsQlFVa3NRMEZCUXl4TFFVRkxMRVZCUVVVN1dVRkRWaXhQUVVGUExFVkJRVVVzUTBGQlF6dFRRVU5ZTzFGQlJVUXNUVUZCVFN4UFFVRlBMRWRCUVVjc1RVRkJRU3hWUVVGVkxHRkJRVllzVlVGQlZTeDFRa0ZCVml4VlFVRlZMRU5CUVVVc1QwRkJUeXh0UTBGQlNTeEpRVUZCTEZOQlFVY3NSVUZCUXl4RFFVRkRMRU5CUVVNc1EwRkJRenRSUVVNNVF5eE5RVUZOTEZkQlFWY3NSMEZCUnl4TlFVRkJMRlZCUVZVc1lVRkJWaXhWUVVGVkxIVkNRVUZXTEZWQlFWVXNRMEZCUlN4WFFVRlhMRzFEUVVGSkxFbEJRVUVzVTBGQlJ5eEZRVU5vUkN4SlFVRkpMRU5CUVVNc1IwRkJSeXhEUVVGRExFZEJRVWNzVFVGQlFTeE5RVUZCTEV0QlFVc3NRMEZCUXl4SlFVRkpMREJEUVVGRkxFZEJRVWNzUTBGQlF5eEhRVUZITEVOQlFVTXNSVUZCUlN4RFFVRkRMRWRCUVVjc1EwRkJReXhOUVVGTkxFTkJRVU1zYlVOQlFVa3NRMEZCUXl4RFFVRkRMRU5CUVVNc1EwRkJReXhGUVVOMFJDeE5RVUZCTEUxQlFVRXNTMEZCU3l4RFFVRkRMRWxCUVVrc01FTkJRVVVzVFVGQlRTeHRRMEZCU1N4RFFVRkRMRU5CUTNoQ0xFTkJRVU03VVVGRFJpeEpRVUZKTEZkQlFWY3NRMEZCUXl4RFFVRkRMRWxCUVVrc1QwRkJUeXhEUVVGRExFTkJRVU1zU1VGQlNTeFhRVUZYTEVOQlFVTXNRMEZCUXl4SlFVRkpMRTlCUVU4c1EwRkJReXhEUVVGRExFVkJRVVU3V1VGRE5VUXNUMEZCVHl4RlFVRkZMRU5CUVVNN1UwRkRXRHRSUVVWRUxFMUJRVTBzVFVGQlRTeEhRVUZuUWl4RlFVRkZMRU5CUVVNN1VVRkRMMElzUzBGQlN5eEpRVUZKTEVOQlFVTXNSMEZCUnl4UFFVRlBMRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zUjBGQlJ5eFhRVUZYTEVOQlFVTXNRMEZCUXl4RlFVRkZMRU5CUVVNc1JVRkJSU3hGUVVGRk8xbEJRemxETEUxQlFVMHNSMEZCUnl4SFFVRmpMRVZCUVVVc1EwRkJRenRaUVVVeFFpeExRVUZMTEVsQlFVa3NRMEZCUXl4SFFVRkhMRTlCUVU4c1EwRkJReXhEUVVGRExFVkJRVVVzUTBGQlF5eEhRVUZITEZkQlFWY3NRMEZCUXl4RFFVRkRMRVZCUVVVc1EwRkJReXhGUVVGRkxFVkJRVVU3WjBKQlF6bERMRTFCUVUwc1VVRkJVU3hIUVVGSExFMUJRVUVzVFVGQlFTeExRVUZMTEVOQlFVTXNTVUZCU1N3d1EwRkJSeXhEUVVGRExFTkJRVU1zTUVOQlFVY3NRMEZCUXl4RFFVRkRMRU5CUVVNN1owSkJRM1JETEVsQlFVa3NVVUZCVVN4TFFVRkxMRk5CUVZNc1NVRkJTU3hSUVVGUkxFdEJRVXNzUTBGQlF5eERRVUZETEVWQlFVVTdiMEpCUXpkRExFZEJRVWNzUTBGQlF5eEpRVUZKTEVOQlFVTXNTMEZCU3l4RFFVRkRMRU5CUVVNN2IwSkJRMmhDTEZOQlFWTTdhVUpCUTFZN1owSkJSVVFzVFVGQlRTeEpRVUZKTEVkQlFVY3NUVUZCUVN4TFFVRkxMRU5CUVVNc1MwRkJTeXd3UTBGQlJ5eFJRVUZSTEVOQlFVTXNRMEZCUXp0blFrRkRja01zU1VGQlNTeERRVUZETEVsQlFVa3NSVUZCUlR0dlFrRkRWQ3hIUVVGSExFTkJRVU1zU1VGQlNTeERRVUZETEV0QlFVc3NRMEZCUXl4RFFVRkRPMjlDUVVOb1FpeFRRVUZUTzJsQ1FVTldPMmRDUVVWRUxFbEJRVWtzVTBGQlV5eEpRVUZKTEVOQlFVTXNTVUZCU3l4RFFVRkRMRk5CUVZNc1EwRkJReXhGUVVGRk8yOUNRVU5zUXl4SFFVRkhMRU5CUVVNc1NVRkJTU3hEUVVGRExFdEJRVXNzUTBGQlF5eERRVUZETzI5Q1FVTm9RaXhUUVVGVE8ybENRVU5XTzJkQ1FVVkVMRWRCUVVjc1EwRkJReXhKUVVGSkxFTkJRVU1zU1VGQlNTeERRVUZETEVOQlFVTTdZVUZEYUVJN1dVRkZSQ3hOUVVGTkxFTkJRVU1zU1VGQlNTeERRVUZETEVkQlFVY3NRMEZCUXl4RFFVRkRPMU5CUTJ4Q08xRkJSVVFzVDBGQlR5eEpRVUZCTEhGRFFVRnJRaXhGUVVGRExFMUJRVTBzUTBGQlF5eERRVUZETzBsQlEzQkRMRU5CUVVNN1NVRkZSRHM3T3pzN096dFBRVTlITzBsQlEwa3NhVUpCUVdsQ0xFTkJRM1JDTEZGQlFXRXNSVUZEWWl4VFFVRnJRanRSUVVWc1FpeEpRVUZKTEZOQlFWTXNSVUZCUlR0WlFVTmlMRTlCUVU4c1NVRkJTU3hEUVVGRExIZENRVUYzUWl4RFFVRkRMRkZCUVZFc1JVRkJSU3hUUVVGVExFTkJRVU1zUTBGQlF6dFRRVU16UkR0UlFVVkVMRTFCUVUwc1RVRkJUU3hIUVVGcFJDeEZRVUZGTEVOQlFVTTdVVUZEYUVVc1MwRkJTeXhOUVVGTkxFdEJRVXNzU1VGQlNTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRTFCUVUwc1JVRkJSVHRaUVVOMlF5eE5RVUZOTEVOQlFVTXNTMEZCU3l4RFFVRkRMRWxCUVVrc1EwRkJReXhIUVVGSExFbEJRVWtzUTBGQlF5eDNRa0ZCZDBJc1EwRkJReXhSUVVGUkxFVkJRVVVzUzBGQlN5eERRVUZETEVsQlFVa3NRMEZCUXl4RFFVRkRPMU5CUXpGRk8xRkJSVVFzVDBGQlR5eE5RVUZOTEVOQlFVTTdTVUZEYUVJc1EwRkJRenRKUVVWUExIZENRVUYzUWl4RFFVTTVRaXhSUVVGaExFVkJRMklzVTBGQmFVSTdPMUZCUldwQ0xFMUJRVTBzV1VGQldTeEhRVUZITEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUXpGQ0xGTkJRVWNzUTBGQlF5eEhRVUZITEVOQlFVTXNVVUZCVVN4RlFVRkZMRU5CUVVNc1IwRkJSeXhKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEZGQlFWRXNRMEZCUXl4RlFVTTFReXhKUVVGSkxFTkJRVU1zUzBGQlN5eERRVU5ZTEVOQlFVTTdVVUZGUml4TlFVRk5MRXRCUVVzc1IwRkJSeXhKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEUxQlFVMHNRMEZCUXl4SlFVRkpMRU5CUVVNc1EwRkJReXhEUVVGRExFVkJRVVVzUlVGQlJTeERRVUZETEVOQlFVTXNRMEZCUXl4SlFVRkpMRXRCUVVzc1UwRkJVeXhEUVVGRExFTkJRVU03VVVGRGNFVXNTVUZCU1N4RFFVRkRMRXRCUVVzc1JVRkJSVHRaUVVOV0xFOUJRVThzU1VGQlNTeERRVUZETzFOQlEySTdVVUZGUkN4TlFVRk5MRkZCUVZFc1IwRkJSeXhOUVVGQkxFMUJRVUVzUzBGQlN5eERRVUZETEVsQlFVa3NNRU5CUVVjc1dVRkJXU3hEUVVGRExFTkJRVU1zUTBGQlF5d3dRMEZCUnl4WlFVRlpMRU5CUVVNc1EwRkJReXhEUVVGRExFTkJRVU03VVVGRGFFVXNTVUZCU1N4UlFVRlJMRXRCUVVzc1UwRkJVeXhKUVVGSkxGRkJRVkVzUzBGQlN5eERRVUZETEVOQlFVTXNSVUZCUlR0WlFVTTNReXhQUVVGUExFbEJRVWtzUTBGQlF6dFRRVU5pTzFGQlJVUXNTVUZCU1N4TFFVRkxMRU5CUVVNc1MwRkJTeXhGUVVGRk8xbEJRMllzVDBGQlR5eE5RVUZCTEV0QlFVc3NRMEZCUXl4TFFVRkxMRU5CUVVNc1VVRkJVU3hEUVVGRExHMURRVUZKTEVsQlFVa3NRMEZCUXp0VFFVTjBRenRSUVVWRUxFOUJRVThzU1VGQlNTeERRVUZETzBsQlEyUXNRMEZCUXp0SlFVVlBMRlZCUVZVc1EwRkJReXhEUVVGTk8xRkJRM1pDTEU5QlFVOHNVMEZCUnl4RFFVRkRMRWRCUVVjc1EwRkJReXhEUVVGRExFTkJRVU1zUTBGQlF6dEpRVU53UWl4RFFVRkRPMGxCUlUwc1NVRkJTU3hEUVVOVUxFOUJRV2xETEVWQlEycERMRTFCUVZjc1JVRkRXQ3hSUVVGaExFVkJRMklzUzBGQllUczdVVUZGWWl4TlFVRk5MR2xDUVVGcFFpeEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hIUVVGSExFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNVMEZCVXl4RFFVRkRPMUZCUTNwRkxFMUJRVTBzVjBGQlZ5eEhRVUZITEVsQlFVRXNVMEZCUnl4RlFVRkRMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVjBGQlZ5eERRVUZETEVOQlFVTTdVVUZGYkVRc2IwSkJRVzlDTzFGQlEzQkNMRWxCUVVrc1YwRkJWeXhIUVVGSExFdEJRVXNzUTBGQlF6dFJRVU40UWl4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVVVGQlVTeEpRVUZKTEZkQlFWY3NSMEZCUnl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGRkJRVkVzUlVGQlJUdFpRVU5vUlN4WFFVRlhMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eFJRVUZSTEVOQlFVTTdVMEZEY2tNN1VVRkRSQ3hKUVVGSkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNVVUZCVVN4SlFVRkpMRmRCUVZjc1IwRkJSeXhKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEZGQlFWRXNSVUZCUlR0WlFVTm9SU3hYUVVGWExFZEJRVWNzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRU5CUVVNN1UwRkRja003VVVGRlJDeHBRMEZCYVVNN1VVRkRha01zU1VGQlNTeGpRVUZqTEVkQlFVY3NTVUZCUVN4VFFVRkhMRVZCUVVNc1VVRkJVU3hEUVVGRExFTkJRVU03VVVGRGJrTXNTVUZCU1N4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFMUJRVTBzU1VGQlNTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMSEZDUVVGeFFpeEZRVUZGTzFsQlF6ZEVMRTFCUVUwc1kwRkJZeXhIUVVGSExFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNVVUZCVVN4SFFVRkhMRmRCUVZjc1EwRkJRenRaUVVNelJDeE5RVUZOTEdkQ1FVRm5RaXhIUVVGSExGTkJRVWNzUTBGQlF5eEhRVUZITEVOQlF6bENMRk5CUVVjc1EwRkJReXhIUVVGSExFTkJRVU1zVFVGQlRTeEZRVUZGTEVOQlFVTXNSMEZCUnl4RFFVRkRMRmRCUVZjc1IwRkJSeXhEUVVGRExFTkJRVU1zUTBGQlF5eEZRVU4wUXl4SlFVRkpMRU5CUVVNc1NVRkJTU3hEUVVOV0xFTkJRVU03V1VGRFJpeE5RVUZOTEZkQlFWY3NSMEZCUnl4SlFVRkJMRk5CUVVjc1JVRkRja0lzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4TlFVRk5MRU5CUVVNc1QwRkJUeXhEUVVGRExFTkJRVU1zUjBGQlJ5eGpRVUZqTEVkQlFVY3NaMEpCUVdkQ0xFTkJRVU1zUTBGQlF5eEZRVU51UlN4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFMUJRVTBzUTBGQlF5eFBRVUZQTEVOQlFVTXNRMEZCUXl4SFFVRkhMR05CUVdNc1IwRkJSeXhuUWtGQlowSXNRMEZCUXl4RFFVRkRMRU5CUTNCRkxFTkJRVU03V1VGRFJpeE5RVUZOTEZkQlFWY3NSMEZCUnl4SlFVRkJMRk5CUVVjc1JVRkRja0lzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4TlFVRk5MRU5CUVVNc1YwRkJWeXhEUVVGRExFTkJRVU1zUjBGQlJ5eGpRVUZqTEVkQlFVY3NaMEpCUVdkQ0xFTkJRVU1zUTBGQlF5eEZRVU4yUlN4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFMUJRVTBzUTBGQlF5eFhRVUZYTEVOQlFVTXNRMEZCUXl4SFFVRkhMR05CUVdNc1IwRkJSeXhuUWtGQlowSXNRMEZCUXl4RFFVRkRMRU5CUTNoRkxFTkJRVU03V1VGRlJpeGpRVUZqTEVkQlFVY3NTVUZCUVN4VFFVRkhMRVZCUTJ4Q0xFbEJRVUVzWVVGQlN5eEZRVUZETEdOQlFXTXNRMEZCUXl4RFFVRkRMRVZCUVVVc1YwRkJWeXhEUVVGRExFTkJRVU1zUlVGQlJTeFhRVUZYTEVOQlFVTXNRMEZCUXl4RFFVRkRMRVZCUTNKRUxFbEJRVUVzWVVGQlN5eEZRVUZETEdOQlFXTXNRMEZCUXl4RFFVRkRMRVZCUVVVc1YwRkJWeXhEUVVGRExFTkJRVU1zUlVGQlJTeFhRVUZYTEVOQlFVTXNRMEZCUXl4RFFVRkRMRU5CUTNSRUxFTkJRVU03VTBGRFNEdFJRVVZFTEUxQlFVMHNhMEpCUVd0Q0xFZEJRVWNzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZEYUVNc1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGRFRDeE5RVUZOTEVWQlEwNHNRMEZCUXl4SFFVRkhMRU5CUVVNc2FVSkJRV2xDTEVkQlFVY3NWMEZCVnl4RFFVRkRMRU5CUTNSRExFVkJRMFFzU1VGQlNTeERRVUZETEVsQlFVa3NRMEZEVml4RFFVRkRPMUZCUTBZc1RVRkJUU3hwUWtGQmFVSXNSMEZCUnl4VFFVRkhMRU5CUVVNc1IwRkJSeXhEUVVNdlFpeFRRVUZITEVOQlFVTXNSMEZCUnl4RFFVRkRMR05CUVdNc1JVRkJSU3hEUVVGRExFZEJRVWNzYVVKQlFXbENMRU5CUVVNc1JVRkRPVU1zU1VGQlNTeERRVUZETEV0QlFVc3NRMEZEV0N4RFFVRkRPMUZCUTBZc1RVRkJUU3haUVVGWkxFZEJRVWNzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZETVVJc1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGRFRDeHBRa0ZCYVVJc1JVRkRha0lzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZEVEN4VFFVRkhMRU5CUVVNc1IwRkJSeXhEUVVGRExHdENRVUZyUWl4RlFVRkZMRWRCUVVjc1EwRkJReXhGUVVOb1F5eEpRVUZKTEVOQlFVTXNTVUZCU1N4RFFVTldMRU5CUTBZc1JVRkRSQ3hYUVVGWExFTkJRMW9zUTBGQlF6dFJRVU5HTEUxQlFVMHNaMEpCUVdkQ0xFZEJRVWNzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZET1VJc1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGRFRDeHBRa0ZCYVVJc1JVRkRha0lzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZEVEN4VFFVRkhMRU5CUVVNc1IwRkJSeXhEUVVGRExHdENRVUZyUWl4RlFVRkZMRWRCUVVjc1EwRkJReXhGUVVOb1F5eEpRVUZKTEVOQlFVTXNTVUZCU1N4RFFVTldMRU5CUTBZc1JVRkRSQ3hYUVVGWExFTkJRMW9zUTBGQlF6dFJRVVZHTEU5QlFVOHNRMEZCUXl4SlFVRkpMRVZCUVVVc1EwRkJRenRSUVVObUxFOUJRVThzUTBGQlF5eExRVUZMTEVOQlFVTXNWMEZCVnl4RlFVRkZMRmRCUVZjc1EwRkJReXhEUVVGRE8xRkJRM2hETEU5QlFVOHNRMEZCUXl4VFFVRlRMRU5CUTJZc1EwRkJReXhqUVVGakxFTkJRVU1zUTBGQlF5eEhRVUZITEUxQlFVMHNRMEZCUXl4RFFVRkRMRWRCUVVjc1EwRkJReXhYUVVGWExFZEJRVWNzUTBGQlF5eERRVUZETEVWQlEyaEVMRU5CUVVNc1kwRkJZeXhEUVVGRExFTkJRVU1zUjBGQlJ5eE5RVUZOTEVOQlFVTXNRMEZCUXl4SFFVRkhMRU5CUVVNc1YwRkJWeXhIUVVGSExFTkJRVU1zUTBGQlF5eERRVU5xUkN4RFFVRkRPMUZCUlVZc1RVRkJRU3hOUVVGQkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVWQlFVTXNVMEZCVXl4dFJFRkRjRUlzVDBGQlR5eEZRVU5RTEVsQlFVa3NSVUZEU2l4TlFVRk5MRVZCUTA0c1kwRkJZeXhGUVVOa0xGZEJRVmNzUTBGRFdpeERRVUZETzFGQlJVWXNaMEpCUVdkQ08xRkJRMmhDTEV0QlFVc3NTVUZCU1N4RFFVRkRMRWRCUVVjc1dVRkJXU3hEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVkQlFVY3NaMEpCUVdkQ0xFTkJRVU1zUTBGQlF5eEZRVUZGTEVOQlFVTXNSVUZCUlN4RlFVRkZPMWxCUTNoRUxFdEJRVXNzU1VGQlNTeERRVUZETEVkQlFVY3NXVUZCV1N4RFFVRkRMRU5CUVVNc1JVRkJSU3hEUVVGRExFZEJRVWNzWjBKQlFXZENMRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zUlVGQlJTeEZRVUZGTzJkQ1FVTjRSQ3hOUVVGTkxHRkJRV0VzUjBGQlJ5eEpRVUZCTEZOQlFVY3NSVUZCUXl4RFFVRkRMRVZCUVVVc1EwRkJReXhEUVVGRExFTkJRVU03WjBKQlEyaERMRTFCUVUwc2NVSkJRWEZDTEVkQlFVY3NVMEZCUnl4RFFVRkRMRWRCUVVjc1EwRkJReXhoUVVGaExFVkJRVVVzYVVKQlFXbENMRU5CUVVNc1EwRkJRenRuUWtGRmVFVXNNa05CUVRKRE8yZENRVU16UXl4TlFVRk5MRk5CUVZNc1IwRkJSeXhKUVVGSkxFTkJRVU1zVlVGQlZTeERRVUZETEdGQlFXRXNRMEZCUXl4RFFVRkRPMmRDUVVOcVJDeEpRVUZKTEVOQlFVTXNTVUZCU1N4RFFVRkRMRmRCUVZjc1EwRkJReXhIUVVGSExFTkJRVU1zVTBGQlV5eERRVUZETEVWQlFVVTdiMEpCUTNCRExFbEJRVWtzUTBGQlF5eFhRVUZYTEVOQlFVTXNSMEZCUnl4RFFVRkRMRk5CUVZNc1JVRkJSU3hKUVVGSkxFTkJRVU1zWVVGQllTeERRVU5vUkN4aFFVRmhMRVZCUTJJc2FVSkJRV2xDTEVOQlEyeENMRU5CUVVNc1EwRkJRenRwUWtGRFNqdG5Ra0ZGUkN4TlFVRk5MRXRCUVVzc1IwRkJSeXhKUVVGSkxFTkJRVU1zVjBGQlZ5eERRVUZETEVkQlFVY3NRMEZCUXl4VFFVRlRMRU5CUVVNc1EwRkJRenRuUWtGRE9VTXNTVUZCU1N4TFFVRkxMRVZCUVVVN2IwSkJRMVFzVDBGQlR5eERRVUZETEZOQlFWTXNRMEZEWml4TFFVRkxMRU5CUVVNc1MwRkJTeXhGUVVOWUxIRkNRVUZ4UWl4RFFVRkRMRU5CUVVNc1JVRkRka0lzY1VKQlFYRkNMRU5CUVVNc1EwRkJReXhEUVVONFFpeERRVUZETzJsQ1FVTklPMkZCUTBZN1UwRkRSanRSUVVWRUxFMUJRVUVzVFVGQlFTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RlFVRkRMRlZCUVZVc2JVUkJRM0pDTEU5QlFVOHNSVUZEVUN4SlFVRkpMRVZCUTBvc1RVRkJUU3hGUVVOT0xHTkJRV01zUlVGRFpDeFhRVUZYTEVOQlExb3NRMEZCUXp0UlFVVkdMSFZDUVVGMVFqdFJRVU4yUWl4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zUzBGQlN5eERRVUZETEdWQlFXVXNSVUZCUlR0WlFVTjBReXhOUVVGTkxGZEJRVmNzUjBGQlJ5eFRRVUZITEVOQlFVTXNSMEZCUnl4RFFVTjZRaXhUUVVGSExFTkJRVU1zUjBGQlJ5eERRVU5NTEdsQ1FVRnBRaXhGUVVOcVFpeFRRVUZITEVOQlFVTXNSMEZCUnl4RFFVTk1MRk5CUVVjc1EwRkJReXhIUVVGSExFTkJRMHdzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZCUXl4clFrRkJhMElzUlVGQlJTeEhRVUZITEVOQlFVTXNSVUZEYUVNc1NVRkJTU3hEUVVGRExFbEJRVWtzUTBGRFZpeEZRVU5FTEVsQlFVRXNVMEZCUnl4RlFVRkRMRU5CUVVNc1EwRkJReXhEUVVOUUxFTkJRMFlzUlVGRFJDeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRk5CUVZNc1EwRkRka0lzUTBGQlF6dFpRVU5HTEUxQlFVMHNaVUZCWlN4SFFVRkhMRk5CUVVjc1EwRkJReXhIUVVGSExFTkJRemRDTEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUTB3c2FVSkJRV2xDTEVWQlEycENMRk5CUVVjc1EwRkJReXhIUVVGSExFTkJRMHdzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZEVEN4VFFVRkhMRU5CUVVNc1IwRkJSeXhEUVVGRExHdENRVUZyUWl4RlFVRkZMRWRCUVVjc1EwRkJReXhGUVVOb1F5eEpRVUZKTEVOQlFVTXNTVUZCU1N4RFFVTldMRVZCUTBRc1NVRkJRU3hUUVVGSExFVkJRVU1zUTBGQlF5eERRVUZETEVOQlExQXNRMEZEUml4RlFVTkVMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVTBGQlV5eERRVU4yUWl4RFFVRkRPMWxCUlVZc1MwRkJTeXhKUVVGSkxFTkJRVU1zUjBGQlJ5eFhRVUZYTEVOQlFVTXNRMEZCUXl4RlFVRkZMRU5CUVVNc1IwRkJSeXhsUVVGbExFTkJRVU1zUTBGQlF5eEZRVUZGTEVOQlFVTXNSVUZCUlN4RlFVRkZPMmRDUVVOMFJDeEpRVUZKTEVOQlFVTXNVVUZCVVN4RFFVTllMRTlCUVU4c1JVRkRVQ3hKUVVGQkxGTkJRVWNzUlVGRFJDeGpRVUZqTEVOQlFVTXNRMEZCUXl4SFFVRkhMRTFCUVUwc1EwRkJReXhEUVVGRExFZEJRVWNzUTBGQlF5eFhRVUZYTEVkQlFVY3NRMEZCUXl4RFFVRkRMRVZCUXk5RExFTkJRVU1zUjBGQlJ5eEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRkZCUVZFc1EwRkRNVUlzUlVGRFJDeEpRVUZCTEZOQlFVY3NSVUZEUkN4alFVRmpMRU5CUVVNc1EwRkJReXhIUVVGSExFMUJRVTBzUTBGQlF5eERRVUZETEVkQlFVY3NRMEZCUXl4WFFVRlhMRWRCUVVjc1EwRkJReXhEUVVGRExFVkJReTlETEVOQlFVTXNSMEZCUnl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGRkJRVkVzUTBGRE1VSXNSVUZEUkN4UFFVRlBMRU5CUVVNc2QwSkJRWGRDTEVWQlEyaERMRTlCUVU4c1EwRkJReXcwUWtGQk5FSXNRMEZEY2tNc1EwRkJRenRoUVVOSU8xbEJRMFFzUzBGQlN5eEpRVUZKTEVOQlFVTXNSMEZCUnl4WFFVRlhMRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zUjBGQlJ5eGxRVUZsTEVOQlFVTXNRMEZCUXl4RlFVRkZMRU5CUVVNc1JVRkJSU3hGUVVGRk8yZENRVU4wUkN4SlFVRkpMRU5CUVVNc1VVRkJVU3hEUVVOWUxFOUJRVThzUlVGRFVDeEpRVUZCTEZOQlFVY3NSVUZEUkN4RFFVRkRMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eFJRVUZSTEVWQlEzcENMR05CUVdNc1EwRkJReXhEUVVGRExFZEJRVWNzVFVGQlRTeERRVUZETEVOQlFVTXNSMEZCUnl4RFFVRkRMRmRCUVZjc1IwRkJSeXhEUVVGRExFTkJRVU1zUTBGRGFFUXNSVUZEUkN4SlFVRkJMRk5CUVVjc1JVRkRSQ3hEUVVGRExFZEJRVWNzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRVZCUTNwQ0xHTkJRV01zUTBGQlF5eERRVUZETEVkQlFVY3NUVUZCVFN4RFFVRkRMRU5CUVVNc1IwRkJSeXhEUVVGRExGZEJRVmNzUjBGQlJ5eERRVUZETEVOQlFVTXNRMEZEYUVRc1JVRkRSQ3hQUVVGUExFTkJRVU1zZDBKQlFYZENMRVZCUTJoRExFOUJRVThzUTBGQlF5dzBRa0ZCTkVJc1EwRkRja01zUTBGQlF6dGhRVU5JTzFOQlEwWTdVVUZGUkN4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zUzBGQlN5eERRVUZETEdkQ1FVRm5RaXhGUVVGRk8xbEJRM1pETEV0QlFVc3NTVUZCU1N4RFFVRkRMRWRCUVVjc1dVRkJXU3hEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVkQlFVY3NaMEpCUVdkQ0xFTkJRVU1zUTBGQlF5eEZRVUZGTEVOQlFVTXNSVUZCUlN4RlFVRkZPMmRDUVVONFJDeEpRVUZKTEVOQlFVTXNVVUZCVVN4RFFVTllMRTlCUVU4c1JVRkRVQ3hKUVVGQkxGTkJRVWNzUlVGRFJDeGpRVUZqTEVOQlFVTXNRMEZCUXl4SFFVRkhMRTFCUVUwc1EwRkJReXhEUVVGRExFZEJRVWNzUTBGQlF5eFhRVUZYTEVkQlFVY3NRMEZCUXl4RFFVRkRMRVZCUXk5RExFTkJRVU1zUjBGQlJ5eHBRa0ZCYVVJc1EwRkRkRUlzUlVGRFJDeEpRVUZCTEZOQlFVY3NSVUZEUkN4alFVRmpMRU5CUVVNc1EwRkJReXhIUVVGSExFMUJRVTBzUTBGQlF5eERRVUZETEVkQlFVY3NRMEZCUXl4WFFVRlhMRWRCUVVjc1EwRkJReXhEUVVGRExFVkJReTlETEVOQlFVTXNSMEZCUnl4cFFrRkJhVUlzUTBGRGRFSXNSVUZEUkN4UFFVRlBMRU5CUVVNc2VVSkJRWGxDTEVWQlEycERMRTlCUVU4c1EwRkJReXcyUWtGQk5rSXNRMEZEZEVNc1EwRkJRenRoUVVOSU8xbEJRMFFzUzBGQlN5eEpRVUZKTEVOQlFVTXNSMEZCUnl4WlFVRlpMRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zUjBGQlJ5eG5Ra0ZCWjBJc1EwRkJReXhEUVVGRExFVkJRVVVzUTBGQlF5eEZRVUZGTEVWQlFVVTdaMEpCUTNoRUxFbEJRVWtzUTBGQlF5eFJRVUZSTEVOQlExZ3NUMEZCVHl4RlFVTlFMRWxCUVVFc1UwRkJSeXhGUVVORUxFTkJRVU1zUjBGQlJ5eHBRa0ZCYVVJc1JVRkRja0lzWTBGQll5eERRVUZETEVOQlFVTXNSMEZCUnl4TlFVRk5MRU5CUVVNc1EwRkJReXhIUVVGSExFTkJRVU1zVjBGQlZ5eEhRVUZITEVOQlFVTXNRMEZCUXl4RFFVTm9SQ3hGUVVORUxFbEJRVUVzVTBGQlJ5eEZRVU5FTEVOQlFVTXNSMEZCUnl4cFFrRkJhVUlzUlVGRGNrSXNZMEZCWXl4RFFVRkRMRU5CUVVNc1IwRkJSeXhOUVVGTkxFTkJRVU1zUTBGQlF5eEhRVUZITEVOQlFVTXNWMEZCVnl4SFFVRkhMRU5CUVVNc1EwRkJReXhEUVVOb1JDeEZRVU5FTEU5QlFVOHNRMEZCUXl4NVFrRkJlVUlzUlVGRGFrTXNUMEZCVHl4RFFVRkRMRFpDUVVFMlFpeERRVU4wUXl4RFFVRkRPMkZCUTBnN1UwRkRSanRSUVVWRUxFbEJRVWtzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4TFFVRkxMRU5CUVVNc1pVRkJaU3hGUVVGRk8xbEJRM1JETEU5QlFVOHNRMEZCUXl4SlFVRkpMRVZCUVVVc1EwRkJRenRaUVVObUxFOUJRVThzUTBGQlF5eFRRVUZUTEVkQlFVY3NUMEZCVHl4RFFVRkRMSGRDUVVGM1FpeERRVUZETzFsQlEzSkVMRTlCUVU4c1EwRkJReXhKUVVGSkxFZEJRVWNzVDBGQlR5eERRVUZETEhOQ1FVRnpRaXhEUVVGRE8xbEJRemxETEU5QlFVOHNRMEZCUXl4WlFVRlpMRWRCUVVjc1VVRkJVU3hEUVVGRE8xbEJRMmhETEU5QlFVOHNRMEZCUXl4VFFVRlRMRWRCUVVjc1VVRkJVU3hEUVVGRE8xbEJSVGRDTEV0QlFVc3NTVUZCU1N4RFFVRkRMRWRCUVVjc1dVRkJXU3hEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVkQlFVY3NaMEpCUVdkQ0xFTkJRVU1zUTBGQlF5eEZRVUZGTEVOQlFVTXNSVUZCUlN4RlFVRkZPMmRDUVVONFJDeExRVUZMTEVsQlFVa3NRMEZCUXl4SFFVRkhMRmxCUVZrc1EwRkJReXhEUVVGRExFVkJRVVVzUTBGQlF5eEhRVUZITEdkQ1FVRm5RaXhEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVWQlFVVXNSVUZCUlR0dlFrRkRlRVFzVDBGQlR5eERRVUZETEZGQlFWRXNRMEZEWkN4SFFVRkhMRU5CUVVNc1MwRkJTeXhEUVVGRExFVkJRVVVzUlVGRFdpeERRVUZETEVkQlFVY3NhVUpCUVdsQ0xFZEJRVWNzYVVKQlFXbENMRWRCUVVjc1EwRkJReXhGUVVNM1F5eERRVUZETEVkQlFVY3NhVUpCUVdsQ0xFZEJRVWNzYVVKQlFXbENMRWRCUVVjc1EwRkJReXhEUVVNNVF5eERRVUZETzJsQ1FVTklPMkZCUTBZN1dVRkZSQ3hQUVVGUExFTkJRVU1zVDBGQlR5eEZRVUZGTEVOQlFVTTdVMEZEYmtJN1VVRkZSQ3hKUVVORkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNTMEZCU3l4RFFVRkRMRlZCUVZVN1dVRkROMElzWjBKQlFXZENMRU5CUVVNc1NVRkJRU3hUUVVGSExFVkJRVU1zUTBGQlF5eEZRVUZGTEVOQlFVTXNRMEZCUXl4RlFVRkZMRmxCUVZrc1JVRkJSU3huUWtGQlowSXNRMEZCUXl4RlFVTXpSRHRaUVVOQkxFbEJRVWtzUTBGQlF5eFRRVUZUTEVOQlExb3NUMEZCVHl4RlFVTlFMRWxCUVVFc1UwRkJSeXhGUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVOQlFVTXNSVUZEVkN4UFFVRlBMRU5CUVVNc2JVSkJRVzFDTEVWQlF6TkNMRTlCUVU4c1EwRkJReXgxUWtGQmRVSXNSVUZETDBJc1QwRkJUeXhEUVVGRExHbENRVUZwUWl4RFFVTXhRaXhEUVVGRE8xTkJRMGc3VVVGRlJDeFBRVUZQTEVOQlFVTXNUMEZCVHl4RlFVRkZMRU5CUVVNN1NVRkRjRUlzUTBGQlF6dEpRVVZQTEdGQlFXRXNRMEZEYmtJc1lVRkJhMElzUlVGRGJFSXNhVUpCUVhsQ096dFJRVVY2UWl4TlFVRk5MRmRCUVZjc1IwRkJSeXhSUVVGUkxFTkJRVU1zWVVGQllTeERRVUZETEZGQlFWRXNRMEZCUXl4RFFVRkRPMUZCUTNKRUxFMUJRVTBzV1VGQldTeEhRVUZITEZkQlFWY3NRMEZCUXl4VlFVRlZMRU5CUVVNc1NVRkJTU3hEUVVGRkxFTkJRVU03VVVGRmJrUXNWMEZCVnl4RFFVRkRMRXRCUVVzc1IwRkJSeXhwUWtGQmFVSXNRMEZCUXp0UlFVTjBReXhYUVVGWExFTkJRVU1zVFVGQlRTeEhRVUZITEdsQ1FVRnBRaXhEUVVGRE8xRkJSWFpETEVsQlFVa3NTMEZCU3l4SFFVRnBRanRaUVVONFFpeGhRVUZoTzFsQlEySXNTMEZCU3l4RlFVRkZMRmRCUVZjN1UwRkRia0lzUTBGQlF6dFJRVVZHTEUxQlFVMHNWMEZCVnl4SFFVRkhMRk5CUVVjc1EwRkJReXhIUVVGSExFTkJRVU1zWVVGQllTeEZRVUZGTEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1UwRkJVeXhEUVVGRExFTkJRVU03VVVGRGJrVXNUVUZCVFN4bFFVRmxMRWRCUVVjc1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGRE4wSXNWMEZCVnl4RlFVTllMRWxCUVVFc1UwRkJSeXhGUVVGRExFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNVMEZCVXl4SFFVRkhMRU5CUVVNc1EwRkJReXhEUVVOb1F5eERRVUZETzFGQlEwWXNUVUZCVFN4aFFVRmhMRWRCUVVjc1RVRkJRU3hOUVVGQkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNUVUZCVFN3d1EwRkJSU3hQUVVGUExHMURRVUZKTEVsQlFVRXNVMEZCUnl4RlFVRkRMRU5CUVVNc1EwRkJReXhEUVVGRE8xRkJSVGRFTEVsQlFVa3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhuUWtGQlowSXNSVUZCUlR0WlFVTnFReXhOUVVGTkxFMUJRVTBzUjBGQlJ5eEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMR2RDUVVGblFpeERRVU14UXl4WlFVRlpMRVZCUTFvc1NVRkJTU3hGUVVOS08yZENRVU5GTEU5QlFVOHNSVUZCUlN4WFFVRlhPMmRDUVVOd1FpeFhRVUZYTEVWQlFVVXNaVUZCWlR0aFFVTTNRaXhGUVVORUxHRkJRV0VzUTBGRFpDeERRVUZETzFsQlJVWXNTVUZCU1N4TFFVRkxMRU5CUVVNc1QwRkJUeXhEUVVGRExFMUJRVTBzUTBGQlF5eEZRVUZGTzJkQ1FVTjZRaXhKUVVGSkxFTkJRVU1zVFVGQlRTeERRVUZETEVOQlFVTXNRMEZCUXl4RlFVRkZPMjlDUVVOa0xFOUJRVThzUzBGQlN5eERRVUZETzJsQ1FVTmtPMkZCUTBZN1UwRkRSanRSUVVWRUxIRkVRVUZ4UkR0UlFVTnlSQ3hMUVVGTExFMUJRVTBzUzBGQlN5eEpRVUZKTEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1RVRkJUU3hGUVVGRk8xbEJRM1pETEZsQlFWa3NRMEZCUXl4SlFVRkpMRVZCUVVVc1EwRkJRenRaUVVOd1FpeFpRVUZaTEVOQlFVTXNWMEZCVnl4SFFVRkhMRTFCUVVFc1MwRkJTeXhEUVVGRExFOUJRVThzYlVOQlFVa3NRMEZCUXl4RFFVRkRPMWxCUlRsRExFMUJRVTBzVTBGQlV5eEhRVUZITEUxQlFVRXNTMEZCU3l4RFFVRkRMRk5CUVZNc2JVTkJRVWtzWVVGQllTeERRVUZETEUxQlFVMHNRMEZCUXp0WlFVVXhSQ3hMUVVGTExFbEJRVWtzUTBGQlF5eEhRVUZITEZkQlFWY3NRMEZCUXl4RFFVRkRMRVZCUVVVc1EwRkJReXhKUVVGSkxHVkJRV1VzUTBGQlF5eERRVUZETEVWQlFVVXNRMEZCUXl4RlFVRkZMRVZCUVVVN1owSkJRM1pFTEV0QlFVc3NTVUZCU1N4RFFVRkRMRWRCUVVjc1YwRkJWeXhEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVsQlFVa3NaVUZCWlN4RFFVRkRMRU5CUVVNc1JVRkJSU3hEUVVGRExFVkJRVVVzUlVGQlJUdHZRa0ZEZGtRc1RVRkJUU3haUVVGWkxFZEJRVWNzU1VGQlFTeFRRVUZITEVWQlFVTXNRMEZCUXl4RlFVRkZMRU5CUVVNc1EwRkJReXhEUVVGRE8yOUNRVVV2UWl4TlFVRkJMRXRCUVVzc1EwRkJReXhoUVVGaExITkVRVU5xUWl4WlFVRlpMRVZCUTFvc1NVRkJTU3hGUVVOS0xFdEJRVXNzUlVGRFRDeGhRVUZoTEVWQlEySXNXVUZCV1N4RFFVTmlMRU5CUVVNN2IwSkJSVVlzVFVGQlRTeG5Ra0ZCWjBJc1IwRkJSeXhUUVVGSExFTkJRVU1zUjBGQlJ5eERRVU01UWl4WlFVRlpMRVZCUTFvc1lVRkJZU3hEUVVOa0xFTkJRVU03YjBKQlJVWXNTVUZCU1N4blFrRkJaMElzUTBGQlF5eERRVUZETEVkQlFVY3NRMEZCUXl4SlFVRkpMR2RDUVVGblFpeERRVUZETEVOQlFVTXNSMEZCUnl4RFFVRkRMRVZCUVVVN2QwSkJRM0JFTEZOQlFWTTdjVUpCUTFZN2IwSkJSVVFzVFVGQlRTeFJRVUZSTEVkQlFVY3NUVUZCUVN4TlFVRkJMRXRCUVVzc1EwRkJReXhKUVVGSkxEQkRRVU4wUWl4blFrRkJaMElzUTBGQlF5eERRVUZETEVOQlFVTXNNRU5CUTI1Q0xHZENRVUZuUWl4RFFVRkRMRU5CUVVNc1EwRkJReXhEUVVGRE8yOUNRVU42UWl4SlFVRkpMRkZCUVZFc1MwRkJTeXhUUVVGVExFbEJRVWtzVVVGQlVTeExRVUZMTEVOQlFVTXNRMEZCUXl4RlFVRkZPM2RDUVVNM1F5eFRRVUZUTzNGQ1FVTldPMjlDUVVWRUxFMUJRVTBzVTBGQlV5eEhRVUZITEUxQlFVRXNUVUZCUVN4TFFVRkxMRU5CUVVNc1MwRkJTeXd3UTBGQlJ5eFJRVUZSTEVOQlFVTXNNRU5CUVVVc1MwRkJTeXhEUVVGRE8yOUNRVU5xUkN4SlFVRkpMRU5CUVVNc1UwRkJVeXhGUVVGRk8zZENRVU5rTEZOQlFWTTdjVUpCUTFZN2IwSkJSVVFzVFVGQlRTeHZRa0ZCYjBJc1IwRkJSeXhUUVVGSExFTkJRVU1zUjBGQlJ5eERRVU5zUXl4VFFVRkhMRU5CUVVNc1IwRkJSeXhEUVVOTUxGbEJRVmtzUlVGRFdpeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRkZCUVZFc1EwRkRkRUlzUlVGRFJDeFRRVUZITEVOQlFVTXNSMEZCUnl4RFFVRkRMR0ZCUVdFc1JVRkJSU3hwUWtGQmFVSXNRMEZCUXl4RFFVTXhReXhEUVVGRE8yOUNRVVZHTEdkQ1FVRm5RanR2UWtGRGFFSXNTVUZCU1N4TFFVRkxMRU5CUVVNc1NVRkJTU3hGUVVGRk8zZENRVU5rTEZsQlFWa3NRMEZCUXl4SlFVRkpMRVZCUVVVc1EwRkJRenQzUWtGRGNFSXNXVUZCV1N4RFFVRkRMRk5CUVZNc1JVRkJSU3hEUVVGRE8zZENRVU42UWl4WlFVRlpMRU5CUVVNc1NVRkJTU3hEUVVObUxHOUNRVUZ2UWl4RFFVRkRMRU5CUVVNc1JVRkRkRUlzYjBKQlFXOUNMRU5CUVVNc1EwRkJReXhGUVVOMFFpeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRkZCUVZFc1JVRkRja0lzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRU5CUTNSQ0xFTkJRVU03ZDBKQlEwWXNXVUZCV1N4RFFVRkRMRWxCUVVrc1JVRkJSU3hEUVVGRE8zRkNRVU55UWp0dlFrRkZSQ3hwUWtGQmFVSTdiMEpCUTJwQ0xFbEJRVWtzZVVKQlFUaENMRU5CUVVNN2IwSkJRMjVETEZGQlFWRXNVMEZCVXl4RlFVRkZPM2RDUVVOcVFpeExRVUZMTEdGQlFXRXNRMEZCUXl4UFFVRlBPelJDUVVONFFpeDVRa0ZCZVVJc1IwRkJSeXhKUVVGQkxGTkJRVWNzUlVGQlF5eHZRa0ZCYjBJc1EwRkJReXhEUVVGRE96UkNRVU4wUkN4TlFVRk5PM2RDUVVWU0xFdEJRVXNzWVVGQllTeERRVUZETEVkQlFVYzdORUpCUTNCQ0xIbENRVUY1UWl4SFFVRkhMRWxCUVVFc1UwRkJSeXhGUVVNM1FpeERRVU5GTEc5Q1FVRnZRaXhEUVVGRExFTkJRVU1zUjBGQlJ5eEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRkZCUVZFc1IwRkJSeXhEUVVGRExFTkJRMjVFTEVkQlFVY3NVMEZCVXl4RFFVRkRMRXRCUVVzc1IwRkJSeXhEUVVGRExFVkJRM1pDTEc5Q1FVRnZRaXhEUVVGRExFTkJRVU1zUTBGRGRrSXNRMEZCUXpzMFFrRkRSaXhOUVVGTk8zZENRVVZTTEV0QlFVc3NZVUZCWVN4RFFVRkRMRkZCUVZFN05FSkJRM3BDTEhsQ1FVRjVRaXhIUVVGSExFbEJRVUVzVTBGQlJ5eEZRVU0zUWl4dlFrRkJiMElzUTBGQlF5eERRVUZETEVkQlFVY3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFZEJRVWNzVTBGQlV5eERRVUZETEV0QlFVc3NSVUZEYUVVc2IwSkJRVzlDTEVOQlFVTXNRMEZCUXl4RFFVTjJRaXhEUVVGRE96UkNRVU5HTEUxQlFVMDdkMEpCUlZJc1MwRkJTeXhoUVVGaExFTkJRVU1zU1VGQlNUczBRa0ZEY2tJc2VVSkJRWGxDTEVkQlFVY3NTVUZCUVN4VFFVRkhMRVZCUXpkQ0xHOUNRVUZ2UWl4RFFVRkRMRU5CUVVNc1JVRkRkRUlzUTBGRFJTeHZRa0ZCYjBJc1EwRkJReXhEUVVGRExFZEJRVWNzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRWRCUVVjc1EwRkJReXhEUVVOdVJDeEhRVUZITEZOQlFWTXNRMEZCUXl4TlFVRk5MRWRCUVVjc1EwRkJReXhEUVVONlFpeERRVUZET3pSQ1FVTkdMRTFCUVUwN2QwSkJSVklzUzBGQlN5eGhRVUZoTEVOQlFVTXNUVUZCVFRzMFFrRkRka0lzZVVKQlFYbENMRWRCUVVjc1NVRkJRU3hUUVVGSExFVkJRemRDTEVOQlEwVXNiMEpCUVc5Q0xFTkJRVU1zUTBGQlF5eEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hIUVVGSExFTkJRVU1zUTBGRGJrUXNSMEZCUnl4VFFVRlRMRU5CUVVNc1MwRkJTeXhIUVVGSExFTkJRVU1zUlVGRGRrSXNRMEZEUlN4dlFrRkJiMElzUTBGQlF5eERRVUZETEVkQlFVY3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFZEJRVWNzUTBGQlF5eERRVU51UkN4SFFVRkhMRk5CUVZNc1EwRkJReXhOUVVGTkxFZEJRVWNzUTBGQlF5eERRVU42UWl4RFFVRkRPelJDUVVOR0xFMUJRVTA3ZDBKQlJWSXNTMEZCU3l4aFFVRmhMRU5CUVVNc1MwRkJTenMwUWtGRGRFSXNlVUpCUVhsQ0xFZEJRVWNzU1VGQlFTeFRRVUZITEVWQlF6ZENMRzlDUVVGdlFpeERRVUZETEVOQlFVTXNSMEZCUnl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGRkJRVkVzUjBGQlJ5eFRRVUZUTEVOQlFVTXNTMEZCU3l4RlFVTm9SU3hEUVVORkxHOUNRVUZ2UWl4RFFVRkRMRU5CUVVNc1IwRkJSeXhKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEZGQlFWRXNSMEZCUnl4RFFVRkRMRU5CUTI1RUxFZEJRVWNzVTBGQlV5eERRVUZETEUxQlFVMHNSMEZCUnl4RFFVRkRMRU5CUTNwQ0xFTkJRVU03TkVKQlEwWXNUVUZCVFR0M1FrRkZVaXhMUVVGTExHRkJRV0VzUTBGQlF5eFZRVUZWT3pSQ1FVTXpRaXg1UWtGQmVVSXNSMEZCUnl4SlFVRkJMRk5CUVVjc1JVRkROMElzYjBKQlFXOUNMRU5CUVVNc1EwRkJReXhGUVVOMFFpeHZRa0ZCYjBJc1EwRkJReXhEUVVGRExFZEJRVWNzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRWRCUVVjc1UwRkJVeXhEUVVGRExFMUJRVTBzUTBGRGJFVXNRMEZCUXpzMFFrRkRSaXhOUVVGTk8zZENRVVZTTEV0QlFVc3NZVUZCWVN4RFFVRkRMRTFCUVUwN05FSkJRM1pDTEhsQ1FVRjVRaXhIUVVGSExFbEJRVUVzVTBGQlJ5eEZRVU0zUWl4RFFVTkZMRzlDUVVGdlFpeERRVUZETEVOQlFVTXNSMEZCUnl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGRkJRVkVzUjBGQlJ5eERRVUZETEVOQlEyNUVMRWRCUVVjc1UwRkJVeXhEUVVGRExFdEJRVXNzUjBGQlJ5eERRVUZETEVWQlEzWkNMRzlDUVVGdlFpeERRVUZETEVOQlFVTXNSMEZCUnl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGRkJRVkVzUjBGQlJ5eFRRVUZUTEVOQlFVTXNUVUZCVFN4RFFVTnNSU3hEUVVGRE96UkNRVU5HTEUxQlFVMDdkMEpCUlZJc1MwRkJTeXhoUVVGaExFTkJRVU1zVjBGQlZ6czBRa0ZETlVJc2VVSkJRWGxDTEVkQlFVY3NTVUZCUVN4VFFVRkhMRVZCUXpkQ0xHOUNRVUZ2UWl4RFFVRkRMRU5CUVVNc1IwRkJSeXhKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEZGQlFWRXNSMEZCUnl4VFFVRlRMRU5CUVVNc1MwRkJTeXhGUVVOb1JTeHZRa0ZCYjBJc1EwRkJReXhEUVVGRExFZEJRVWNzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRWRCUVVjc1UwRkJVeXhEUVVGRExFMUJRVTBzUTBGRGJFVXNRMEZCUXpzMFFrRkRSaXhOUVVGTk8zRkNRVU5VTzI5Q1FVVkVMRmxCUVZrc1EwRkJReXhUUVVGVExFTkJRM0JDTEZOQlFWTXNSVUZEVkN4NVFrRkJlVUlzUTBGQlF5eERRVUZETEVWQlF6TkNMSGxDUVVGNVFpeERRVUZETEVOQlFVTXNRMEZETlVJc1EwRkJRenR2UWtGRlJpeEpRVUZKTEV0QlFVc3NRMEZCUXl4SlFVRkpMRVZCUVVVN2QwSkJRMlFzV1VGQldTeERRVUZETEU5QlFVOHNSVUZCUlN4RFFVRkRPM0ZDUVVONFFqdHZRa0ZGUkN4TlFVRkJMRXRCUVVzc1EwRkJReXhqUVVGakxITkVRVU5zUWl4WFFVRlhMRVZCUTFnc1dVRkJXU3hGUVVOYUxFbEJRVWtzUlVGRFNpeExRVUZMTEVWQlEwd3NZVUZCWVN4RlFVTmlMRmxCUVZrc1EwRkRZaXhEUVVGRE8ybENRVU5JTzJGQlEwWTdXVUZGUkN4WlFVRlpMRU5CUVVNc1QwRkJUeXhGUVVGRkxFTkJRVU03VTBGRGVFSTdVVUZGUkN4TlFVRkJMRTFCUVVFc1NVRkJTU3hEUVVGRExFOUJRVThzUlVGQlF5eHBRa0ZCYVVJc2JVUkJRelZDTEZkQlFWY3NSVUZEV0N4WlFVRlpMRVZCUTFvc1NVRkJTU3hGUVVOS08xbEJRMFVzVDBGQlR5eEZRVUZGTEZkQlFWYzdXVUZEY0VJc1YwRkJWeXhGUVVGRkxHVkJRV1U3VTBGRE4wSXNSVUZEUkN4aFFVRmhMRU5CUTJRc1EwRkJRenRSUVVWR0xFOUJRVThzUzBGQlN5eERRVUZETzBsQlEyWXNRMEZCUXp0SlFVVlBMRkZCUVZFc1EwRkRaQ3hQUVVGcFF5eEZRVU5xUXl4TFFVRlZMRVZCUTFZc1IwRkJVU3hGUVVOU0xFMUJRV01zUlVGRFpDeFRRVUZwUWp0UlFVVnFRaXhQUVVGUExFTkJRVU1zU1VGQlNTeEZRVUZGTEVOQlFVTTdVVUZGWml4UFFVRlBMRU5CUVVNc1UwRkJVeXhIUVVGSExGTkJRVk1zUTBGQlF6dFJRVU01UWl4UFFVRlBMRU5CUVVNc1YwRkJWeXhIUVVGSExFMUJRVTBzUTBGQlF6dFJRVVUzUWl4UFFVRlBMRU5CUVVNc1UwRkJVeXhGUVVGRkxFTkJRVU03VVVGRGNFSXNUMEZCVHl4RFFVRkRMRTFCUVUwc1EwRkJReXhMUVVGTExFTkJRVU1zUTBGQlF5eEZRVUZGTEV0QlFVc3NRMEZCUXl4RFFVRkRMRU5CUVVNc1EwRkJRenRSUVVOcVF5eFBRVUZQTEVOQlFVTXNUVUZCVFN4RFFVRkRMRWRCUVVjc1EwRkJReXhEUVVGRExFVkJRVVVzUjBGQlJ5eERRVUZETEVOQlFVTXNRMEZCUXl4RFFVRkRPMUZCUXpkQ0xFOUJRVThzUTBGQlF5eE5RVUZOTEVWQlFVVXNRMEZCUXp0UlFVVnFRaXhQUVVGUExFTkJRVU1zVDBGQlR5eEZRVUZGTEVOQlFVTTdTVUZEY0VJc1EwRkJRenRKUVVWUExGTkJRVk1zUTBGRFppeFBRVUZwUXl4RlFVTnFReXhSUVVGaExFVkJRMklzVFVGQll5eEZRVU5rTEZOQlFXbENMRVZCUTJwQ0xFbEJRVms3VVVGRldpeFBRVUZQTEVOQlFVTXNTVUZCU1N4RlFVRkZMRU5CUVVNN1VVRkZaaXhQUVVGUExFTkJRVU1zVTBGQlV5eEhRVUZITEZOQlFWTXNRMEZCUXp0UlFVVTVRaXhOUVVGTkxGRkJRVkVzUjBGQlJ5eEpRVUZKTEVOQlFVTXNTVUZCU1N4RFFVRkRMRWxCUVVrc1IwRkJSeXhEUVVGRExFTkJRVU1zUTBGQlF6dFJRVU55UXl4UFFVRlBMRU5CUVVNc1YwRkJWeXhIUVVGSExFMUJRVTBzUTBGQlF6dFJRVU0zUWl4UFFVRlBMRU5CUVVNc1UwRkJVeXhGUVVGRkxFTkJRVU03VVVGRGNFSXNUMEZCVHl4RFFVRkRMRTFCUVUwc1EwRkJReXhSUVVGUkxFTkJRVU1zUTBGQlF5eEhRVUZITEZGQlFWRXNSVUZCUlN4UlFVRlJMRU5CUVVNc1EwRkJReXhIUVVGSExGRkJRVkVzUTBGQlF5eERRVUZETzFGQlF6ZEVMRTlCUVU4c1EwRkJReXhOUVVGTkxFTkJRVU1zVVVGQlVTeERRVUZETEVOQlFVTXNSMEZCUnl4UlFVRlJMRVZCUVVVc1VVRkJVU3hEUVVGRExFTkJRVU1zUjBGQlJ5eFJRVUZSTEVOQlFVTXNRMEZCUXp0UlFVTTNSQ3hQUVVGUExFTkJRVU1zVFVGQlRTeERRVUZETEZGQlFWRXNRMEZCUXl4RFFVRkRMRWRCUVVjc1VVRkJVU3hGUVVGRkxGRkJRVkVzUTBGQlF5eERRVUZETEVkQlFVY3NVVUZCVVN4RFFVRkRMRU5CUVVNN1VVRkROMFFzVDBGQlR5eERRVUZETEUxQlFVMHNRMEZCUXl4UlFVRlJMRU5CUVVNc1EwRkJReXhIUVVGSExGRkJRVkVzUlVGQlJTeFJRVUZSTEVOQlFVTXNRMEZCUXl4SFFVRkhMRkZCUVZFc1EwRkJReXhEUVVGRE8xRkJRemRFTEU5QlFVOHNRMEZCUXl4TlFVRk5MRVZCUVVVc1EwRkJRenRSUVVWcVFpeFBRVUZQTEVOQlFVTXNUMEZCVHl4RlFVRkZMRU5CUVVNN1NVRkRjRUlzUTBGQlF6czdRVUZ3Y1VKSUxEQkNRWEZ4UWtNN1FVRndjVUo1UWl4MVFrRkJaU3hIUVVGdFFqdEpRVU40UkN4eFFrRkJjVUlzUlVGQlJTeEpRVUZKTzBsQlF6TkNMRkZCUVZFc1JVRkJSU3hGUVVGRk8wbEJRMW9zVFVGQlRTeEZRVUZGTzFGQlEwNDdXVUZEUlN4SlFVRkpMRVZCUVVVc1UwRkJVenRUUVVOb1FqdExRVU5HTzBsQlEwUXNVMEZCVXl4RlFVRkZMRU5CUVVNN1NVRkRXaXhYUVVGWExFVkJRVVVzUTBGQlF6dEpRVU5rTEd0Q1FVRnJRaXhGUVVGRkxFVkJRVVU3UTBGRGRrSXNRMEZCUXp0QlFVVnpRaXd5UWtGQmJVSXNSMEZCUnl4TlFVRk5MRU5CUVVNN1FVRkROMElzSzBKQlFYVkNMRWRCUVVjc1EwRkJReXhEUVVGRE8wRkJRelZDTEhsQ1FVRnBRaXhIUVVGSExFVkJRVVVzUTBGQlF6dEJRVVYyUWl4cFEwRkJlVUlzUjBGQlJ5eFJRVUZSTEVOQlFVTTdRVUZEY2tNc2NVTkJRVFpDTEVkQlFVY3NRMEZCUXl4RFFVRkRPMEZCUld4RExHZERRVUYzUWl4SFFVRkhMRTlCUVU4c1EwRkJRenRCUVVOdVF5dzRRa0ZCYzBJc1IwRkJSeXhuUWtGQlowSXNRMEZCUXp0QlFVVXhReXhuUTBGQmQwSXNSMEZCUnl4UlFVRlJMRU5CUVVNN1FVRkRjRU1zYjBOQlFUUkNMRWRCUVVjc1EwRkJReXhEUVVGRE8wRkJPRzlDTTBRN096czdPMGRCUzBjN1FVRkRTU3hMUVVGTExGVkJRVlVzT0VKQlFUaENMRU5CUTJ4RUxFOUJTMFVzUlVGRFJpeEpRVXRETEVWQlEwUXNUMEZGUlR0SlFVVkdMRTFCUVUwc2JVSkJRVzFDTEVkQlFVY3NRMEZCUXl4SlFVRlpMRVZCUjJoRExFVkJRVVU3TzFGQlExUXNUVUZCVFN4TFFVRkxMRWRCUVVjc1RVRkJRU3hQUVVGUExFTkJRVU1zU1VGQlNTeERRVUZETERCRFFVRkZMRTlCUVU4c1EwRkJRenRSUVVOeVF5eEpRVUZKTEVOQlFVTXNTMEZCU3l4RlFVRkZPMWxCUTFZc1RVRkJUU3hKUVVGSkxFdEJRVXNzUTBGQlF5eFZRVUZWTEVsQlFVa3NZVUZCWVN4RFFVRkRMRU5CUVVNN1UwRkRPVU03VVVGRlJDeFBRVUZQTEV0QlFVc3NRMEZCUXp0SlFVTm1MRU5CUVVNc1EwRkJRenRKUVVWR0xFMUJRVTBzVFVGQlRTeEhRVUZSTEVsQlFVa3NRMEZCUXp0SlFVTjZRaXhKUVVGSkxFMUJRVTBzUTBGQlF5eE5RVUZOTEVWQlFVVTdVVUZEYWtJc1MwRkJTeXhOUVVGTkxFTkJRVU1zUTBGQlF5eEZRVUZGTEV0QlFVc3NRMEZCUXl4SlFVRkpMRTFCUVUwc1EwRkJReXhOUVVGTkxFTkJRVU1zVDBGQlR5eEZRVUZGTEVWQlFVVTdXVUZEYUVRc05rUkJRVFpFTzFsQlF6ZEVMRWxCUVVrc1MwRkJTeXhEUVVGRExFdEJRVXNzUlVGQlJUdG5Ra0ZEWml4TFFVRkxMRTFCUVUwc1EwRkJReXhEUVVGRExFVkJRVVVzU1VGQlNTeERRVUZETEVsQlFVa3NTMEZCU3l4RFFVRkRMRXRCUVVzc1EwRkJReXhQUVVGUExFVkJRVVVzUlVGQlJUdHZRa0ZETjBNc1RVRkJUU3hEUVVGRExFMUJRVTBzUTBGQlF5eERRVUZETEVOQlFVTXNRMEZCUXl4TFFVRkxMRU5CUVVNc1EwRkJReXhEUVVGRExFTkJRVU1zUzBGQlN5eEhRVUZITEcxQ1FVRnRRaXhEUVVGRExFbEJRVWtzUTBGQlF5eFRRVUZUTEVOQlFVTXNRMEZCUXp0dlFrRkRkRVVzVDBGQlR5eE5RVUZOTEVOQlFVTXNUVUZCVFN4RFFVRkRMRU5CUVVNc1EwRkJReXhEUVVGRExFdEJRVXNzUTBGQlF5eERRVUZETEVOQlFVTXNRMEZCUXl4VFFVRlRMRU5CUVVNN2FVSkJRelZETzJGQlEwWTdXVUZGUkN4M1FrRkJkMEk3V1VGRGVFSXNTVUZCU1N4RFFVRkJMRTlCUVU4c1lVRkJVQ3hQUVVGUExIVkNRVUZRTEU5QlFVOHNRMEZCUlN4alFVRmpMRXRCUVVrc1MwRkJTeXhEUVVGRExFbEJRVWtzU1VGQlNTeExRVUZMTEVOQlFVTXNTMEZCU3l4RlFVRkZPMmRDUVVONFJDeE5RVUZOTEVOQlFVTXNUVUZCVFN4RFFVRkRMRU5CUVVNc1EwRkJReXhEUVVGRExFbEJRVWtzUjBGQlJ5eEpRVUZCTEdGQlFVc3NSVUZCUXl4SlFVRkJMR2RDUVVGTkxFVkJRVU1zUzBGQlN5eERRVUZETEVsQlFVa3NRMEZCUXl4RlFVRkZMRXRCUVVzc1EwRkJReXhMUVVGTExFTkJRVU1zUTBGQlF6dG5Ra0ZETDBRc1QwRkJUeXhOUVVGTkxFTkJRVU1zVFVGQlRTeERRVUZETEVOQlFVTXNRMEZCUXl4RFFVRkRMRXRCUVVzc1EwRkJRenRoUVVNdlFqdFRRVU5HTzB0QlEwWTdTVUZGUkN4aFFVRmhPMGxCUTJJc1NVRkJTU3hEUVVGRExFOUJRVThzUjBGQlJ5eE5RVUV5UWl4RFFVRkRPMEZCUXpkRExFTkJRVU03UVVGc1JFUXNkMFZCYTBSREluMD0iLCIvLyBUaGUgbW9kdWxlIGNhY2hlXG52YXIgX193ZWJwYWNrX21vZHVsZV9jYWNoZV9fID0ge307XG5cbi8vIFRoZSByZXF1aXJlIGZ1bmN0aW9uXG5mdW5jdGlvbiBfX3dlYnBhY2tfcmVxdWlyZV9fKG1vZHVsZUlkKSB7XG5cdC8vIENoZWNrIGlmIG1vZHVsZSBpcyBpbiBjYWNoZVxuXHR2YXIgY2FjaGVkTW9kdWxlID0gX193ZWJwYWNrX21vZHVsZV9jYWNoZV9fW21vZHVsZUlkXTtcblx0aWYgKGNhY2hlZE1vZHVsZSAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0cmV0dXJuIGNhY2hlZE1vZHVsZS5leHBvcnRzO1xuXHR9XG5cdC8vIENyZWF0ZSBhIG5ldyBtb2R1bGUgKGFuZCBwdXQgaXQgaW50byB0aGUgY2FjaGUpXG5cdHZhciBtb2R1bGUgPSBfX3dlYnBhY2tfbW9kdWxlX2NhY2hlX19bbW9kdWxlSWRdID0ge1xuXHRcdC8vIG5vIG1vZHVsZS5pZCBuZWVkZWRcblx0XHQvLyBubyBtb2R1bGUubG9hZGVkIG5lZWRlZFxuXHRcdGV4cG9ydHM6IHt9XG5cdH07XG5cblx0Ly8gRXhlY3V0ZSB0aGUgbW9kdWxlIGZ1bmN0aW9uXG5cdF9fd2VicGFja19tb2R1bGVzX19bbW9kdWxlSWRdLmNhbGwobW9kdWxlLmV4cG9ydHMsIG1vZHVsZSwgbW9kdWxlLmV4cG9ydHMsIF9fd2VicGFja19yZXF1aXJlX18pO1xuXG5cdC8vIFJldHVybiB0aGUgZXhwb3J0cyBvZiB0aGUgbW9kdWxlXG5cdHJldHVybiBtb2R1bGUuZXhwb3J0cztcbn1cblxuIiwiIiwiLy8gc3RhcnR1cFxuLy8gTG9hZCBlbnRyeSBtb2R1bGUgYW5kIHJldHVybiBleHBvcnRzXG4vLyBUaGlzIGVudHJ5IG1vZHVsZSBpcyByZWZlcmVuY2VkIGJ5IG90aGVyIG1vZHVsZXMgc28gaXQgY2FuJ3QgYmUgaW5saW5lZFxudmFyIF9fd2VicGFja19leHBvcnRzX18gPSBfX3dlYnBhY2tfcmVxdWlyZV9fKFwiLi9pbmRleC50c1wiKTtcbiIsIiJdLCJuYW1lcyI6W10sInNvdXJjZVJvb3QiOiIifQ==