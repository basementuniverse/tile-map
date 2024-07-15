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

/***/ "./node_modules/lru_map/dist/lru.js":
/*!******************************************!*\
  !*** ./node_modules/lru_map/dist/lru.js ***!
  \******************************************/
/***/ (function(__unused_webpack_module, exports) {

!function(g,c){ true?c(exports):0}(this,function(g){const c=Symbol("newer"),e=Symbol("older");class n{constructor(a,b){typeof a!=="number"&&(b=a,a=0),this.size=0,this.limit=a,this.oldest=this.newest=void 0,this._keymap=new Map(),b&&(this.assign(b),a<1&&(this.limit=this.size))}_markEntryAsUsed(a){if(a===this.newest)return;a[c]&&(a===this.oldest&&(this.oldest=a[c]),a[c][e]=a[e]),a[e]&&(a[e][c]=a[c]),a[c]=void 0,a[e]=this.newest,this.newest&&(this.newest[c]=a),this.newest=a}assign(a){let b,d=this.limit||Number.MAX_VALUE;this._keymap.clear();let m=a[Symbol.iterator]();for(let h=m.next();!h.done;h=m.next()){let f=new l(h.value[0],h.value[1]);this._keymap.set(f.key,f),b?(b[c]=f,f[e]=b):this.oldest=f,b=f;if(d--==0)throw new Error("overflow")}this.newest=b,this.size=this._keymap.size}get(a){var b=this._keymap.get(a);return b?(this._markEntryAsUsed(b),b.value):void 0}set(a,b){var d=this._keymap.get(a);return d?(d.value=b,this._markEntryAsUsed(d),this):(this._keymap.set(a,d=new l(a,b)),this.newest?(this.newest[c]=d,d[e]=this.newest):this.oldest=d,this.newest=d,++this.size,this.size>this.limit&&this.shift(),this)}shift(){var a=this.oldest;if(a)return this.oldest[c]?(this.oldest=this.oldest[c],this.oldest[e]=void 0):(this.oldest=void 0,this.newest=void 0),a[c]=a[e]=void 0,this._keymap.delete(a.key),--this.size,[a.key,a.value]}find(a){let b=this._keymap.get(a);return b?b.value:void 0}has(a){return this._keymap.has(a)}delete(a){var b=this._keymap.get(a);return b?(this._keymap.delete(b.key),b[c]&&b[e]?(b[e][c]=b[c],b[c][e]=b[e]):b[c]?(b[c][e]=void 0,this.oldest=b[c]):b[e]?(b[e][c]=void 0,this.newest=b[e]):this.oldest=this.newest=void 0,this.size--,b.value):void 0}clear(){this.oldest=this.newest=void 0,this.size=0,this._keymap.clear()}keys(){return new j(this.oldest)}values(){return new k(this.oldest)}entries(){return this}[Symbol.iterator](){return new i(this.oldest)}forEach(a,b){typeof b!=="object"&&(b=this);let d=this.oldest;for(;d;)a.call(b,d.value,d.key,this),d=d[c]}toJSON(){for(var a=new Array(this.size),b=0,d=this.oldest;d;)a[b++]={key:d.key,value:d.value},d=d[c];return a}toString(){for(var a="",b=this.oldest;b;)a+=String(b.key)+":"+b.value,b=b[c],b&&(a+=" < ");return a}}g.LRUMap=n;function l(a,b){this.key=a,this.value=b,this[c]=void 0,this[e]=void 0}function i(a){this.entry=a}i.prototype[Symbol.iterator]=function(){return this},i.prototype.next=function(){let a=this.entry;return a?(this.entry=a[c],{done:!1,value:[a.key,a.value]}):{done:!0,value:void 0}};function j(a){this.entry=a}j.prototype[Symbol.iterator]=function(){return this},j.prototype.next=function(){let a=this.entry;return a?(this.entry=a[c],{done:!1,value:a.key}):{done:!0,value:void 0}};function k(a){this.entry=a}k.prototype[Symbol.iterator]=function(){return this},k.prototype.next=function(){let a=this.entry;return a?(this.entry=a[c],{done:!1,value:a.value}):{done:!0,value:void 0}}});
//# sourceMappingURL=lru.js.map


/***/ }),

/***/ "./index.ts":
/*!******************!*\
  !*** ./index.ts ***!
  \******************/
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";

Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.tileMapOptionsContentProcessor = exports.TileMap = exports.TileAlignment = void 0;
const vec_1 = __webpack_require__(/*! @basementuniverse/vec */ "./node_modules/@basementuniverse/vec/vec.js");
const lru_map_1 = __webpack_require__(/*! lru_map */ "./node_modules/lru_map/dist/lru.js");
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
function clamp(a, min = 0, max = 1) {
    return a < min ? min : (a > max ? max : a);
}
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
     * Get a minimal set of rectangles which cover the tiles in a given layer
     *
     * @param layerName The name of the layer to get rectangles for
     * @param fieldName We will check the truthyness of this field in the
     * tile definition
     * @param tileBounds Optional bounds to check
     */
    getLayerRectangles(layerName, fieldName, tileBounds) {
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
            actualPosition = (0, vec_1.vec)(clamp(actualPosition.x, minPosition.x, maxPosition.x), clamp(actualPosition.y, minPosition.y, maxPosition.y));
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
// TODO
// https://www.npmjs.com/package/rectangle-decomposition
// https://www.npmjs.com/package/fast-rle
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
async function tileMapOptionsContentProcessor(content, data) {
    //
}
exports.tileMapOptionsContentProcessor = tileMapOptionsContentProcessor;
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiaW5kZXguanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyIuLi9pbmRleC50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7QUFBQSwrQ0FBNEM7QUFDNUMscUNBQWlDO0FBNlFqQyxJQUFZLGFBVVg7QUFWRCxXQUFZLGFBQWE7SUFDdkIsdURBQVcsQ0FBQTtJQUNYLCtDQUFHLENBQUE7SUFDSCx5REFBUSxDQUFBO0lBQ1IsaURBQUksQ0FBQTtJQUNKLHFEQUFNLENBQUE7SUFDTixtREFBSyxDQUFBO0lBQ0wsNkRBQVUsQ0FBQTtJQUNWLHFEQUFNLENBQUE7SUFDTiwrREFBVyxDQUFBO0FBQ2IsQ0FBQyxFQVZXLGFBQWEsR0FBYixxQkFBYSxLQUFiLHFCQUFhLFFBVXhCO0FBY0QsU0FBUyxLQUFLLENBQUMsQ0FBUyxFQUFFLE1BQWMsQ0FBQyxFQUFFLE1BQWMsQ0FBQztJQUN4RCxPQUFPLENBQUMsR0FBRyxHQUFHLENBQUMsQ0FBQyxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLEdBQUcsR0FBRyxDQUFDLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQzdDLENBQUM7QUFFRCxTQUFTLGdCQUFnQixDQUN2QixLQUFVLEVBQ1YsT0FBWSxFQUNaLFdBQWdCO0lBRWhCLE9BQU8sQ0FDTCxLQUFLLENBQUMsQ0FBQyxJQUFJLE9BQU8sQ0FBQyxDQUFDO1FBQ3BCLEtBQUssQ0FBQyxDQUFDLElBQUksT0FBTyxDQUFDLENBQUM7UUFDcEIsS0FBSyxDQUFDLENBQUMsR0FBRyxXQUFXLENBQUMsQ0FBQztRQUN2QixLQUFLLENBQUMsQ0FBQyxHQUFHLFdBQVcsQ0FBQyxDQUFDLENBQ3hCLENBQUM7QUFDSixDQUFDO0FBRUQsTUFBYSxPQUFPO0lBc0NsQixZQUFtQixPQUFvQztRQUNyRCxNQUFNLGFBQWEsR0FBRyxNQUFNLENBQUMsTUFBTSxDQUNqQyxFQUFFLEVBQ0YsT0FBTyxDQUFDLGVBQWUsRUFDdkIsT0FB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

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
//# sourceMappingURL=data:application/json;charset=utf-8;base64,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Q0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHNDQUFzQztBQUN0QywyQ0FBMkM7Ozs7OztVQzFVM0M7VUFDQTs7VUFFQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTs7VUFFQTtVQUNBOztVQUVBO1VBQ0E7VUFDQTs7OztVRXRCQTtVQUNBO1VBQ0E7VUFDQSIsInNvdXJjZXMiOlsid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3RpbGUtbWFwL3dlYnBhY2svdW5pdmVyc2FsTW9kdWxlRGVmaW5pdGlvbiIsIndlYnBhY2s6Ly9AYmFzZW1lbnR1bml2ZXJzZS90aWxlLW1hcC8uL25vZGVfbW9kdWxlcy9AYmFzZW1lbnR1bml2ZXJzZS91dGlscy91dGlscy5qcyIsIndlYnBhY2s6Ly9AYmFzZW1lbnR1bml2ZXJzZS90aWxlLW1hcC8uL25vZGVfbW9kdWxlcy9AYmFzZW1lbnR1bml2ZXJzZS92ZWMvdmVjLmpzIiwid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3RpbGUtbWFwLy4vbm9kZV9tb2R1bGVzL2xydV9tYXAvZGlzdC9scnUuanMiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2UvdGlsZS1tYXAvLi9pbmRleC50cyIsIndlYnBhY2s6Ly9AYmFzZW1lbnR1bml2ZXJzZS90aWxlLW1hcC93ZWJwYWNrL2Jvb3RzdHJhcCIsIndlYnBhY2s6Ly9AYmFzZW1lbnR1bml2ZXJzZS90aWxlLW1hcC93ZWJwYWNrL2JlZm9yZS1zdGFydHVwIiwid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3RpbGUtbWFwL3dlYnBhY2svc3RhcnR1cCIsIndlYnBhY2s6Ly9AYmFzZW1lbnR1bml2ZXJzZS90aWxlLW1hcC93ZWJwYWNrL2FmdGVyLXN0YXJ0dXAiXSwic291cmNlc0NvbnRlbnQiOlsiKGZ1bmN0aW9uIHdlYnBhY2tVbml2ZXJzYWxNb2R1bGVEZWZpbml0aW9uKHJvb3QsIGZhY3RvcnkpIHtcblx0aWYodHlwZW9mIGV4cG9ydHMgPT09ICdvYmplY3QnICYmIHR5cGVvZiBtb2R1bGUgPT09ICdvYmplY3QnKVxuXHRcdG1vZHVsZS5leHBvcnRzID0gZmFjdG9yeSgpO1xuXHRlbHNlIGlmKHR5cGVvZiBkZWZpbmUgPT09ICdmdW5jdGlvbicgJiYgZGVmaW5lLmFtZClcblx0XHRkZWZpbmUoW10sIGZhY3RvcnkpO1xuXHRlbHNlIHtcblx0XHR2YXIgYSA9IGZhY3RvcnkoKTtcblx0XHRmb3IodmFyIGkgaW4gYSkgKHR5cGVvZiBleHBvcnRzID09PSAnb2JqZWN0JyA/IGV4cG9ydHMgOiByb290KVtpXSA9IGFbaV07XG5cdH1cbn0pKHNlbGYsICgpID0+IHtcbnJldHVybiAiLCIvKipcbiAqIEBvdmVydmlldyBBIGxpYnJhcnkgb2YgdXNlZnVsIGZ1bmN0aW9uc1xuICogQGF1dGhvciBHb3Jkb24gTGFycmlnYW5cbiAqL1xuXG4vKipcbiAqIENoZWNrIGlmIHR3byBudW1iZXJzIGFyZSBhcHByb3hpbWF0ZWx5IGVxdWFsXG4gKiBAcGFyYW0ge251bWJlcn0gYSBOdW1iZXIgYVxuICogQHBhcmFtIHtudW1iZXJ9IGIgTnVtYmVyIGJcbiAqIEBwYXJhbSB7bnVtYmVyfSBbcD1OdW1iZXIuRVBTSUxPTl0gVGhlIHByZWNpc2lvbiB2YWx1ZVxuICogQHJldHVybiB7Ym9vbGVhbn0gVHJ1ZSBpZiBudW1iZXJzIGEgYW5kIGIgYXJlIGFwcHJveGltYXRlbHkgZXF1YWxcbiAqL1xuY29uc3QgZmxvYXRFcXVhbHMgPSAoYSwgYiwgcCA9IE51bWJlci5FUFNJTE9OKSA9PiBNYXRoLmFicyhhIC0gYikgPCBwO1xuXG4vKipcbiAqIENsYW1wIGEgbnVtYmVyIGJldHdlZW4gbWluIGFuZCBtYXhcbiAqIEBwYXJhbSB7bnVtYmVyfSBhIFRoZSBudW1iZXIgdG8gY2xhbXBcbiAqIEBwYXJhbSB7bnVtYmVyfSBbbWluPTBdIFRoZSBtaW5pbXVtIHZhbHVlXG4gKiBAcGFyYW0ge251bWJlcn0gW21heD0xXSBUaGUgbWF4aW11bSB2YWx1ZVxuICogQHJldHVybiB7bnVtYmVyfSBBIGNsYW1wZWQgbnVtYmVyXG4gKi9cbmNvbnN0IGNsYW1wID0gKGEsIG1pbiA9IDAsIG1heCA9IDEpID0+IGEgPCBtaW4gPyBtaW4gOiAoYSA+IG1heCA/IG1heCA6IGEpO1xuXG4vKipcbiAqIEdldCB0aGUgZnJhY3Rpb25hbCBwYXJ0IG9mIGEgbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gYSBUaGUgbnVtYmVyIGZyb20gd2hpY2ggdG8gZ2V0IHRoZSBmcmFjdGlvbmFsIHBhcnRcbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIGZyYWN0aW9uYWwgcGFydCBvZiB0aGUgbnVtYmVyXG4gKi9cbmNvbnN0IGZyYWMgPSBhID0+IGEgPj0gMCA/IGEgLSBNYXRoLmZsb29yKGEpIDogYSAtIE1hdGguY2VpbChhKTtcblxuLyoqXG4gKiBSb3VuZCBuIHRvIGQgZGVjaW1hbCBwbGFjZXNcbiAqIEBwYXJhbSB7bnVtYmVyfSBuIFRoZSBudW1iZXIgdG8gcm91bmRcbiAqIEBwYXJhbSB7bnVtYmVyfSBbZD0wXSBUaGUgbnVtYmVyIG9mIGRlY2ltYWwgcGxhY2VzIHRvIHJvdW5kIHRvXG4gKiBAcmV0dXJuIHtudW1iZXJ9IEEgcm91bmRlZCBudW1iZXJcbiAqL1xuY29uc3Qgcm91bmQgPSAobiwgZCA9IDApID0+IHtcbiAgY29uc3QgcCA9IE1hdGgucG93KDEwLCBkKTtcbiAgcmV0dXJuIE1hdGgucm91bmQobiAqIHAgKyBOdW1iZXIuRVBTSUxPTikgLyBwO1xufVxuXG4vKipcbiAqIERvIGEgbGluZWFyIGludGVycG9sYXRpb24gYmV0d2VlbiBhIGFuZCBiXG4gKiBAcGFyYW0ge251bWJlcn0gYSBUaGUgbWluaW11bSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBiIFRoZSBtYXhpbXVtIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIGludGVycG9sYXRpb24gdmFsdWUsIHNob3VsZCBiZSBpbiB0aGUgaW50ZXJ2YWwgWzAsIDFdXG4gKiBAcmV0dXJuIHtudW1iZXJ9IEFuIGludGVycG9sYXRlZCB2YWx1ZSBpbiB0aGUgaW50ZXJ2YWwgW2EsIGJdXG4gKi9cbmNvbnN0IGxlcnAgPSAoYSwgYiwgaSkgPT4gYSArIChiIC0gYSkgKiBpO1xuXG4vKipcbiAqIEdldCB0aGUgcG9zaXRpb24gb2YgaSBiZXR3ZWVuIGEgYW5kIGJcbiAqIEBwYXJhbSB7bnVtYmVyfSBhIFRoZSBtaW5pbXVtIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGIgVGhlIG1heGltdW0gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgaW50ZXJwb2xhdGVkIHZhbHVlIGluIHRoZSBpbnRlcnZhbCBbYSwgYl1cbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIHBvc2l0aW9uIG9mIGkgYmV0d2VlbiBhIGFuZCBiXG4gKi9cbmNvbnN0IHVubGVycCA9IChhLCBiLCBpKSA9PiAoaSAtIGEpIC8gKGIgLSBhKTtcblxuLyoqXG4gKiBEbyBhIGJpbGluZWFyIGludGVycG9sYXRpb25cbiAqIEBwYXJhbSB7bnVtYmVyfSBjMDAgVG9wLWxlZnQgdmFsdWVcbiAqIEBwYXJhbSB7bnVtYmVyfSBjMTAgVG9wLXJpZ2h0IHZhbHVlXG4gKiBAcGFyYW0ge251bWJlcn0gYzAxIEJvdHRvbS1sZWZ0IHZhbHVlXG4gKiBAcGFyYW0ge251bWJlcn0gYzExIEJvdHRvbS1yaWdodCB2YWx1ZVxuICogQHBhcmFtIHtudW1iZXJ9IGl4IEludGVycG9sYXRpb24gdmFsdWUgYWxvbmcgeFxuICogQHBhcmFtIHtudW1iZXJ9IGl5IEludGVycG9sYXRpb24gdmFsdWUgYWxvbmcgeVxuICogQHJldHVybiB7bnVtYmVyfSBBIGJpbGluZWFyIGludGVycG9sYXRlZCB2YWx1ZVxuICovXG5jb25zdCBibGVycCA9IChjMDAsIGMxMCwgYzAxLCBjMTEsIGl4LCBpeSkgPT4gbGVycChsZXJwKGMwMCwgYzEwLCBpeCksIGxlcnAoYzAxLCBjMTEsIGl4KSwgaXkpO1xuXG4vKipcbiAqIFJlLW1hcCBhIG51bWJlciBpIGZyb20gcmFuZ2UgYTEuLi5hMiB0byBiMS4uLmIyXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgbnVtYmVyIHRvIHJlLW1hcFxuICogQHBhcmFtIHtudW1iZXJ9IGExXG4gKiBAcGFyYW0ge251bWJlcn0gYTJcbiAqIEBwYXJhbSB7bnVtYmVyfSBiMVxuICogQHBhcmFtIHtudW1iZXJ9IGIyXG4gKiBAcmV0dXJuIHtudW1iZXJ9XG4gKi9cbmNvbnN0IHJlbWFwID0gKGksIGExLCBhMiwgYjEsIGIyKSA9PiBiMSArIChpIC0gYTEpICogKGIyIC0gYjEpIC8gKGEyIC0gYTEpO1xuXG4vKipcbiAqIERvIGEgc21vb3RoIGludGVycG9sYXRpb24gYmV0d2VlbiBhIGFuZCBiXG4gKiBAcGFyYW0ge251bWJlcn0gYSBUaGUgbWluaW11bSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBiIFRoZSBtYXhpbXVtIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIGludGVycG9sYXRpb24gdmFsdWVcbiAqIEByZXR1cm4ge251bWJlcn0gQW4gaW50ZXJwb2xhdGVkIHZhbHVlIGluIHRoZSBpbnRlcnZhbCBbYSwgYl1cbiAqL1xuY29uc3Qgc21vb3Roc3RlcCA9IChhLCBiLCBpKSA9PiBsZXJwKGEsIGIsIDMgKiBNYXRoLnBvdyhpLCAyKSAtIDIgKiBNYXRoLnBvdyhpLCAzKSk7XG5cbi8qKlxuICogR2V0IGFuIGFuZ2xlIGluIHJhZGlhbnNcbiAqIEBwYXJhbSB7bnVtYmVyfSBkZWdyZWVzIFRoZSBhbmdsZSBpbiBkZWdyZWVzXG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSBhbmdsZSBpbiByYWRpYW5zXG4gKi9cbmNvbnN0IHJhZGlhbnMgPSBkZWdyZWVzID0+IChNYXRoLlBJIC8gMTgwKSAqIGRlZ3JlZXM7XG5cbi8qKlxuICogR2V0IGFuIGFuZ2xlIGluIGRlZ3JlZXNcbiAqIEBwYXJhbSB7bnVtYmVyfSByYWRpYW5zIFRoZSBhbmdsZSBpbiByYWRpYW5zXG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSBhbmdsZSBpbiBkZWdyZWVzXG4gKi9cbmNvbnN0IGRlZ3JlZXMgPSByYWRpYW5zID0+ICgxODAgLyBNYXRoLlBJKSAqIHJhZGlhbnM7XG5cbi8qKlxuICogR2V0IGEgcmFuZG9tIGZsb2F0IGluIHRoZSBpbnRlcnZhbCBbbWluLCBtYXgpXG4gKiBAcGFyYW0ge251bWJlcn0gbWluIEluY2x1c2l2ZSBtaW5cbiAqIEBwYXJhbSB7bnVtYmVyfSBtYXggRXhjbHVzaXZlIG1heFxuICogQHJldHVybiB7bnVtYmVyfSBBIHJhbmRvbSBmbG9hdCBpbiB0aGUgaW50ZXJ2YWwgW21pbiwgbWF4KVxuICovXG5jb25zdCByYW5kb21CZXR3ZWVuID0gKG1pbiwgbWF4KSA9PiBNYXRoLnJhbmRvbSgpICogKG1heCAtIG1pbikgKyBtaW47XG5cbi8qKlxuICogR2V0IGEgcmFuZG9tIGludGVnZXIgaW4gdGhlIGludGVydmFsIFttaW4sIG1heF1cbiAqIEBwYXJhbSB7bnVtYmVyfSBtaW4gSW5jbHVzaXZlIG1pblxuICogQHBhcmFtIHtudW1iZXJ9IG1heCBJbmNsdXNpdmUgbWF4XG4gKiBAcmV0dXJuIHtudW1iZXJ9IEEgcmFuZG9tIGludGVnZXIgaW4gdGhlIGludGVydmFsIFttaW4sIG1heF1cbiAqL1xuY29uc3QgcmFuZG9tSW50QmV0d2VlbiA9IChtaW4sIG1heCkgPT4gTWF0aC5mbG9vcihNYXRoLnJhbmRvbSgpICogKG1heCAtIG1pbiArIDEpKSArIG1pbjtcblxuLyoqXG4gKiBHZXQgYSBub3JtYWxseS1kaXN0cmlidXRlZCByYW5kb20gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gW211PTAuNV0gVGhlIG1lYW4gdmFsdWVcbiAqIEBwYXJhbSB7bnVtYmVyfSBbc2lnbWE9MC41XSBUaGUgc3RhbmRhcmQgZGV2aWF0aW9uXG4gKiBAcGFyYW0ge251bWJlcn0gW3NhbXBsZXM9Ml0gVGhlIG51bWJlciBvZiBzYW1wbGVzXG4gKiBAcmV0dXJuIHtudW1iZXJ9IEEgbm9ybWFsbHktZGlzdHJpYnV0ZWQgcmFuZG9tIG51bWJlclxuICovXG5jb25zdCBjbHRSYW5kb20gPSAobXUgPSAwLjUsIHNpZ21hID0gMC41LCBzYW1wbGVzID0gMikgPT4ge1xuICBsZXQgdG90YWwgPSAwO1xuICBmb3IgKGxldCBpID0gc2FtcGxlczsgaS0tOykge1xuICAgIHRvdGFsICs9IE1hdGgucmFuZG9tKCk7XG4gIH1cbiAgcmV0dXJuIG11ICsgKHRvdGFsIC0gc2FtcGxlcyAvIDIpIC8gKHNhbXBsZXMgLyAyKSAqIHNpZ21hO1xufTtcblxuLyoqXG4gKiBHZXQgYSBub3JtYWxseS1kaXN0cmlidXRlZCByYW5kb20gaW50ZWdlciBpbiB0aGUgaW50ZXJ2YWwgW21pbiwgbWF4XVxuICogQHBhcmFtIHtudW1iZXJ9IG1pbiBJbmNsdXNpdmUgbWluXG4gKiBAcGFyYW0ge251bWJlcn0gbWF4IEluY2x1c2l2ZSBtYXhcbiAqIEByZXR1cm4ge251bWJlcn0gQSBub3JtYWxseS1kaXN0cmlidXRlZCByYW5kb20gaW50ZWdlclxuICovXG5jb25zdCBjbHRSYW5kb21JbnQgPSAobWluLCBtYXgpID0+IE1hdGguZmxvb3IobWluICsgY2x0UmFuZG9tKDAuNSwgMC41LCAyKSAqIChtYXggKyAxIC0gbWluKSk7XG5cbi8qKlxuICogUmV0dXJuIGEgd2VpZ2h0ZWQgcmFuZG9tIGludGVnZXJcbiAqIEBwYXJhbSB7QXJyYXk8bnVtYmVyPn0gdyBBbiBhcnJheSBvZiB3ZWlnaHRzXG4gKiBAcmV0dXJuIHtudW1iZXJ9IEFuIGluZGV4IGZyb20gd1xuICovXG5jb25zdCB3ZWlnaHRlZFJhbmRvbSA9IHcgPT4ge1xuICBsZXQgdG90YWwgPSB3LnJlZHVjZSgoYSwgaSkgPT4gYSArIGksIDApLCBuID0gMDtcbiAgY29uc3QgciA9IE1hdGgucmFuZG9tKCkgKiB0b3RhbDtcbiAgd2hpbGUgKHRvdGFsID4gcikge1xuICAgIHRvdGFsIC09IHdbbisrXTtcbiAgfVxuICByZXR1cm4gbiAtIDE7XG59O1xuXG4vKipcbiAqIEFuIGludGVycG9sYXRpb24gZnVuY3Rpb25cbiAqIEBjYWxsYmFjayBJbnRlcnBvbGF0aW9uRnVuY3Rpb25cbiAqIEBwYXJhbSB7bnVtYmVyfSBhIFRoZSBtaW5pbXVtIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGIgVGhlIG1heGltdW0gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgaW50ZXJwb2xhdGlvbiB2YWx1ZSwgc2hvdWxkIGJlIGluIHRoZSBpbnRlcnZhbCBbMCwgMV1cbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIGludGVycG9sYXRlZCB2YWx1ZSBpbiB0aGUgaW50ZXJ2YWwgW2EsIGJdXG4gKi9cblxuLyoqXG4gKiBSZXR1cm4gYW4gaW50ZXJwb2xhdGVkIHZhbHVlIGZyb20gYW4gYXJyYXlcbiAqIEBwYXJhbSB7QXJyYXk8bnVtYmVyPn0gYSBBbiBhcnJheSBvZiB2YWx1ZXMgaW50ZXJwb2xhdGVcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIEEgbnVtYmVyIGluIHRoZSBpbnRlcnZhbCBbMCwgMV1cbiAqIEBwYXJhbSB7SW50ZXJwb2xhdGlvbkZ1bmN0aW9ufSBbZj1NYXRoLmxlcnBdIFRoZSBpbnRlcnBvbGF0aW9uIGZ1bmN0aW9uIHRvIHVzZVxuICogQHJldHVybiB7bnVtYmVyfSBBbiBpbnRlcnBvbGF0ZWQgdmFsdWUgaW4gdGhlIGludGVydmFsIFttaW4oYSksIG1heChhKV1cbiAqL1xuY29uc3QgbGVycEFycmF5ID0gKGEsIGksIGYgPSBsZXJwKSA9PiB7XG4gIGNvbnN0IHMgPSBpICogKGEubGVuZ3RoIC0gMSk7XG4gIGNvbnN0IHAgPSBjbGFtcChNYXRoLnRydW5jKHMpLCAwLCBhLmxlbmd0aCAtIDEpO1xuICByZXR1cm4gZihhW3BdIHx8IDAsIGFbcCArIDFdIHx8IDAsIGZyYWMocykpO1xufTtcblxuLyoqXG4gKiBHZXQgdGhlIGRvdCBwcm9kdWN0IG9mIHR3byB2ZWN0b3JzXG4gKiBAcGFyYW0ge0FycmF5PG51bWJlcj59IGEgVmVjdG9yIGFcbiAqIEBwYXJhbSB7QXJyYXk8bnVtYmVyPn0gYiBWZWN0b3IgYlxuICogQHJldHVybiB7bnVtYmVyfSBhIOKImSBiXG4gKi9cbmNvbnN0IGRvdCA9IChhLCBiKSA9PiBhLnJlZHVjZSgobiwgdiwgaSkgPT4gbiArIHYgKiBiW2ldLCAwKTtcblxuLyoqXG4gKiBHZXQgdGhlIGZhY3RvcmlhbCBvZiBhIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGFcbiAqIEByZXR1cm4ge251bWJlcn0gYSFcbiAqL1xuY29uc3QgZmFjdG9yaWFsID0gYSA9PiB7XG4gIGxldCByZXN1bHQgPSAxO1xuICBmb3IgKGxldCBpID0gMjsgaSA8PSBhOyBpKyspIHtcbiAgICByZXN1bHQgKj0gaTtcbiAgfVxuICByZXR1cm4gcmVzdWx0O1xufTtcblxuLyoqXG4gKiBHZXQgdGhlIG51bWJlciBvZiBwZXJtdXRhdGlvbnMgb2YgciBlbGVtZW50cyBmcm9tIGEgc2V0IG9mIG4gZWxlbWVudHNcbiAqIEBwYXJhbSB7bnVtYmVyfSBuXG4gKiBAcGFyYW0ge251bWJlcn0gclxuICogQHJldHVybiB7bnVtYmVyfSBuUHJcbiAqL1xuY29uc3QgbnByID0gKG4sIHIpID0+IGZhY3RvcmlhbChuKSAvIGZhY3RvcmlhbChuIC0gcik7XG5cbi8qKlxuICogR2V0IHRoZSBudW1iZXIgb2YgY29tYmluYXRpb25zIG9mIHIgZWxlbWVudHMgZnJvbSBhIHNldCBvZiBuIGVsZW1lbnRzXG4gKiBAcGFyYW0ge251bWJlcn0gblxuICogQHBhcmFtIHtudW1iZXJ9IHJcbiAqIEByZXR1cm4ge251bWJlcn0gbkNyXG4gKi9cbmNvbnN0IG5jciA9IChuLCByKSA9PiBmYWN0b3JpYWwobikgLyAoZmFjdG9yaWFsKHIpICogZmFjdG9yaWFsKG4gLSByKSk7XG5cbi8qKlxuICogR2VuZXJhdGUgYWxsIGNvbWJpbmF0aW9ucyBvZiByIGVsZW1lbnRzIGZyb20gYW4gYXJyYXlcbiAqXG4gKiBAZXhhbXBsZVxuICogYGBganNcbiAqIGNvbWJpbmF0aW9ucyhbMSwgMiwgM10sIDIpO1xuICogYGBgXG4gKlxuICogT3V0cHV0OlxuICogYGBganNvblxuICogW1xuICogICBbMSwgMl0sXG4gKiAgIFsxLCAzXSxcbiAqICAgWzIsIDNdXG4gKiBdXG4gKiBgYGBcbiAqIEBwYXJhbSB7QXJyYXk8Kj59IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSByIFRoZSBudW1iZXIgb2YgZWxlbWVudHMgdG8gY2hvb3NlIGluIGVhY2ggY29tYmluYXRpb25cbiAqIEByZXR1cm4ge0FycmF5PEFycmF5PCo+Pn0gQW4gYXJyYXkgb2YgY29tYmluYXRpb24gYXJyYXlzXG4gKi9cbmNvbnN0IGNvbWJpbmF0aW9ucyA9IChhLCByKSA9PiB7XG4gIGlmIChyID09PSAxKSB7XG4gICAgcmV0dXJuIGEubWFwKGl0ZW0gPT4gW2l0ZW1dKTtcbiAgfVxuXG4gIHJldHVybiBhLnJlZHVjZShcbiAgICAoYWNjLCBpdGVtLCBpKSA9PiBbXG4gICAgICAuLi5hY2MsXG4gICAgICAuLi5jb21iaW5hdGlvbnMoYS5zbGljZShpICsgMSksIHIgLSAxKS5tYXAoYyA9PiBbaXRlbSwgLi4uY10pLFxuICAgIF0sXG4gICAgW11cbiAgKTtcbn07XG5cbi8qKlxuICogR2V0IGEgY2FydGVzaWFuIHByb2R1Y3Qgb2YgYXJyYXlzXG4gKlxuICogQGV4YW1wbGVcbiAqIGBgYGpzXG4gKiBjYXJ0ZXNpYW4oWzEsIDIsIDNdLCBbJ2EnLCAnYiddKTtcbiAqIGBgYFxuICpcbiAqIE91dHB1dDpcbiAqIGBgYGpzb25cbiAqIFtcbiAqICAgWzEsIFwiYVwiXSxcbiAqICAgWzEsIFwiYlwiXSxcbiAqICAgWzIsIFwiYVwiXSxcbiAqICAgWzIsIFwiYlwiXSxcbiAqICAgWzMsIFwiYVwiXSxcbiAqICAgWzMsIFwiYlwiXVxuICogXVxuICogYGBgXG4gKi9cbmNvbnN0IGNhcnRlc2lhbiA9ICguLi5hcnIpID0+XG4gIGFyci5yZWR1Y2UoXG4gICAgKGEsIGIpID0+IGEuZmxhdE1hcChjID0+IGIubWFwKGQgPT4gWy4uLmMsIGRdKSksXG4gICAgW1tdXVxuICApO1xuXG4vKipcbiAqIEEgZnVuY3Rpb24gZm9yIGdlbmVyYXRpbmcgYXJyYXkgdmFsdWVzXG4gKiBAY2FsbGJhY2sgVGltZXNGdW5jdGlvblxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIGFycmF5IGluZGV4XG4gKiBAcmV0dXJuIHsqfSBUaGUgYXJyYXkgdmFsdWVcbiAqL1xuXG4vKipcbiAqIFJldHVybiBhIG5ldyBhcnJheSB3aXRoIGxlbmd0aCBuIGJ5IGNhbGxpbmcgZnVuY3Rpb24gZihpKSBvbiBlYWNoIGVsZW1lbnRcbiAqIEBwYXJhbSB7VGltZXNGdW5jdGlvbn0gZlxuICogQHBhcmFtIHtudW1iZXJ9IG4gVGhlIHNpemUgb2YgdGhlIGFycmF5XG4gKiBAcmV0dXJuIHtBcnJheTwqPn1cbiAqL1xuY29uc3QgdGltZXMgPSAoZiwgbikgPT4gQXJyYXkobikuZmlsbCgwKS5tYXAoKF8sIGkpID0+IGYoaSkpO1xuXG4vKipcbiAqIFJldHVybiBhbiBhcnJheSBjb250YWluaW5nIG51bWJlcnMgMC0+KG4gLSAxKVxuICogQHBhcmFtIHtudW1iZXJ9IG4gVGhlIHNpemUgb2YgdGhlIGFycmF5XG4gKiBAcmV0dXJuIHtBcnJheTxudW1iZXI+fSBBbiBhcnJheSBvZiBpbnRlZ2VycyAwLT4obiAtIDEpXG4gKi9cbmNvbnN0IHJhbmdlID0gbiA9PiB0aW1lcyhpID0+IGksIG4pO1xuXG4vKipcbiAqIFppcCAyIGFycmF5cyB0b2dldGhlciwgaS5lLiAoWzEsIDIsIDNdLCBbYSwgYiwgY10pID0+IFtbMSwgYV0sIFsyLCBiXSwgWzMsIGNdXVxuICogQHBhcmFtIHtBcnJheTwqPn0gYVxuICogQHBhcmFtIHtBcnJheTwqPn0gYlxuICogQHJldHVybiB7QXJyYXk8QXJyYXk8Kj4+fVxuICovXG5jb25zdCB6aXAgPSAoYSwgYikgPT4gYS5tYXAoKGssIGkpID0+IFtrLCBiW2ldXSk7XG5cbi8qKlxuICogUmV0dXJuIGFycmF5W2ldIHdpdGggcG9zaXRpdmUgYW5kIG5lZ2F0aXZlIHdyYXBwaW5nXG4gKiBAcGFyYW0ge0FycmF5PCo+fSBhXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgcG9zaXRpdmVseS9uZWdhdGl2ZWx5IHdyYXBwZWQgYXJyYXkgaW5kZXhcbiAqIEByZXR1cm4geyp9IEFuIGVsZW1lbnQgZnJvbSB0aGUgYXJyYXlcbiAqL1xuY29uc3QgYXQgPSAoYSwgaSkgPT4gYVtpIDwgMCA/IGEubGVuZ3RoIC0gKE1hdGguYWJzKGkgKyAxKSAlIGEubGVuZ3RoKSAtIDEgOiBpICUgYS5sZW5ndGhdO1xuXG4vKipcbiAqIFJldHVybiB0aGUgbGFzdCBlbGVtZW50IG9mIGFuIGFycmF5IHdpdGhvdXQgcmVtb3ZpbmcgaXRcbiAqIEBwYXJhbSB7QXJyYXk8Kj59IGFcbiAqIEByZXR1cm4geyp9IFRoZSBsYXN0IGVsZW1lbnQgZnJvbSB0aGUgYXJyYXlcbiAqL1xuY29uc3QgcGVlayA9IChhKSA9PiB7XG4gIGlmICghYS5sZW5ndGgpIHtcbiAgICByZXR1cm4gdW5kZWZpbmVkO1xuICB9XG5cbiAgcmV0dXJuIGFbYS5sZW5ndGggLSAxXTtcbn07XG5cbi8qKlxuICogQ2hvcCBhbiBhcnJheSBpbnRvIGNodW5rcyBvZiBzaXplIG5cbiAqIEBwYXJhbSB7QXJyYXk8Kj59IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBuIFRoZSBjaHVuayBzaXplXG4gKiBAcmV0dXJuIHtBcnJheTxBcnJheTwqPj59IEFuIGFycmF5IG9mIGFycmF5IGNodW5rc1xuICovXG5jb25zdCBjaHVuayA9IChhLCBuKSA9PiB0aW1lcyhpID0+IGEuc2xpY2UoaSAqIG4sIGkgKiBuICsgbiksIE1hdGguY2VpbChhLmxlbmd0aCAvIG4pKTtcblxuLyoqXG4gKiBSYW5kb21seSBzaHVmZmxlIGEgc2hhbGxvdyBjb3B5IG9mIGFuIGFycmF5XG4gKiBAcGFyYW0ge0FycmF5PCo+fSBhXG4gKiBAcmV0dXJuIHtBcnJheTwqPn0gVGhlIHNodWZmbGVkIGFycmF5XG4gKi9cbmNvbnN0IHNodWZmbGUgPSBhID0+IGEuc2xpY2UoKS5zb3J0KCgpID0+IE1hdGgucmFuZG9tKCkgLSAwLjUpO1xuXG4vKipcbiAqIEZsYXR0ZW4gYW4gb2JqZWN0XG4gKiBAcGFyYW0ge29iamVjdH0gb1xuICogQHBhcmFtIHtzdHJpbmd9IGNvbmNhdGVuYXRvciBUaGUgc3RyaW5nIHRvIHVzZSBmb3IgY29uY2F0ZW5hdGluZyBrZXlzXG4gKiBAcmV0dXJuIHtvYmplY3R9IEEgZmxhdHRlbmVkIG9iamVjdFxuICovXG5jb25zdCBmbGF0ID0gKG8sIGNvbmNhdGVuYXRvciA9ICcuJykgPT4ge1xuICByZXR1cm4gT2JqZWN0LmtleXMobykucmVkdWNlKChhY2MsIGtleSkgPT4ge1xuICAgIGlmIChvW2tleV0gaW5zdGFuY2VvZiBEYXRlKSB7XG4gICAgICByZXR1cm4ge1xuICAgICAgICAuLi5hY2MsXG4gICAgICAgIFtrZXldOiBvW2tleV0udG9JU09TdHJpbmcoKSxcbiAgICAgIH07XG4gICAgfVxuXG4gICAgaWYgKHR5cGVvZiBvW2tleV0gIT09ICdvYmplY3QnIHx8ICFvW2tleV0pIHtcbiAgICAgIHJldHVybiB7XG4gICAgICAgIC4uLmFjYyxcbiAgICAgICAgW2tleV06IG9ba2V5XSxcbiAgICAgIH07XG4gICAgfVxuICAgIGNvbnN0IGZsYXR0ZW5lZCA9IGZsYXQob1trZXldLCBjb25jYXRlbmF0b3IpO1xuXG4gICAgcmV0dXJuIHtcbiAgICAgIC4uLmFjYyxcbiAgICAgIC4uLk9iamVjdC5rZXlzKGZsYXR0ZW5lZCkucmVkdWNlKFxuICAgICAgICAoY2hpbGRBY2MsIGNoaWxkS2V5KSA9PiAoe1xuICAgICAgICAgIC4uLmNoaWxkQWNjLFxuICAgICAgICAgIFtgJHtrZXl9JHtjb25jYXRlbmF0b3J9JHtjaGlsZEtleX1gXTogZmxhdHRlbmVkW2NoaWxkS2V5XSxcbiAgICAgICAgfSksXG4gICAgICAgIHt9XG4gICAgICApLFxuICAgIH07XG4gIH0sIHt9KTtcbn07XG5cbi8qKlxuICogVW5mbGF0dGVuIGFuIG9iamVjdFxuICogQHBhcmFtIHtvYmplY3R9IG9cbiAqIEBwYXJhbSB7c3RyaW5nfSBjb25jYXRlbmF0b3IgVGhlIHN0cmluZyB0byBjaGVjayBmb3IgaW4gY29uY2F0ZW5hdGVkIGtleXNcbiAqIEByZXR1cm4ge29iamVjdH0gQW4gdW4tZmxhdHRlbmVkIG9iamVjdFxuICovXG5jb25zdCB1bmZsYXQgPSAobywgY29uY2F0ZW5hdG9yID0gJy4nKSA9PiB7XG4gIGxldCByZXN1bHQgPSB7fSwgdGVtcCwgc3Vic3RyaW5ncywgcHJvcGVydHksIGk7XG5cbiAgZm9yIChwcm9wZXJ0eSBpbiBvKSB7XG4gICAgc3Vic3RyaW5ncyA9IHByb3BlcnR5LnNwbGl0KGNvbmNhdGVuYXRvcik7XG4gICAgdGVtcCA9IHJlc3VsdDtcbiAgICBmb3IgKGkgPSAwOyBpIDwgc3Vic3RyaW5ncy5sZW5ndGggLSAxOyBpKyspIHtcbiAgICAgIGlmICghKHN1YnN0cmluZ3NbaV0gaW4gdGVtcCkpIHtcbiAgICAgICAgaWYgKGlzRmluaXRlKHN1YnN0cmluZ3NbaSArIDFdKSkge1xuICAgICAgICAgIHRlbXBbc3Vic3RyaW5nc1tpXV0gPSBbXTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICB0ZW1wW3N1YnN0cmluZ3NbaV1dID0ge307XG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHRlbXAgPSB0ZW1wW3N1YnN0cmluZ3NbaV1dO1xuICAgIH1cbiAgICB0ZW1wW3N1YnN0cmluZ3Nbc3Vic3RyaW5ncy5sZW5ndGggLSAxXV0gPSBvW3Byb3BlcnR5XTtcbiAgfVxuXG4gIHJldHVybiByZXN1bHQ7XG59O1xuXG4vKipcbiAqIEEgc3BsaXQgcHJlZGljYXRlXG4gKiBAY2FsbGJhY2sgU3BsaXRQcmVkaWNhdGVcbiAqIEBwYXJhbSB7YW55fSB2YWx1ZSBUaGUgY3VycmVudCB2YWx1ZVxuICogQHJldHVybiB7Ym9vbGVhbn0gVHJ1ZSBpZiB0aGUgYXJyYXkgc2hvdWxkIHNwbGl0IGF0IHRoaXMgaW5kZXhcbiAqL1xuXG4vKipcbiAqIFNwbGl0IGFuIGFycmF5IGludG8gc3ViLWFycmF5cyBiYXNlZCBvbiBhIHByZWRpY2F0ZVxuICogQHBhcmFtIHtBcnJheTwqPn0gYXJyYXlcbiAqIEBwYXJhbSB7U3BsaXRQcmVkaWNhdGV9IHByZWRpY2F0ZVxuICogQHJldHVybiB7QXJyYXk8QXJyYXk8Kj4+fSBBbiBhcnJheSBvZiBhcnJheXNcbiAqL1xuY29uc3Qgc3BsaXQgPSAoYXJyYXksIHByZWRpY2F0ZSkgPT4ge1xuICBjb25zdCByZXN1bHQgPSBbXTtcbiAgbGV0IGN1cnJlbnQgPSBbXTtcbiAgZm9yIChjb25zdCB2YWx1ZSBvZiBhcnJheSkge1xuICAgIGlmIChwcmVkaWNhdGUodmFsdWUpKSB7XG4gICAgICBpZiAoY3VycmVudC5sZW5ndGgpIHtcbiAgICAgICAgcmVzdWx0LnB1c2goY3VycmVudCk7XG4gICAgICB9XG4gICAgICBjdXJyZW50ID0gW3ZhbHVlXTtcbiAgICB9IGVsc2Uge1xuICAgICAgY3VycmVudC5wdXNoKHZhbHVlKTtcbiAgICB9XG4gIH1cbiAgcmVzdWx0LnB1c2goY3VycmVudCk7XG5cbiAgcmV0dXJuIHJlc3VsdDtcbn07XG5cbi8qKlxuICogUGx1Y2sga2V5cyBmcm9tIGFuIG9iamVjdFxuICogQHBhcmFtIHtvYmplY3R9IG9cbiAqIEBwYXJhbSB7Li4uc3RyaW5nfSBrZXlzIFRoZSBrZXlzIHRvIHBsdWNrIGZyb20gdGhlIG9iamVjdFxuICogQHJldHVybiB7b2JqZWN0fSBBbiBvYmplY3QgY29udGFpbmluZyB0aGUgcGx1Y2tlZCBrZXlzXG4gKi9cbmNvbnN0IHBsdWNrID0gKG8sIC4uLmtleXMpID0+IHtcbiAgcmV0dXJuIGtleXMucmVkdWNlKFxuICAgIChyZXN1bHQsIGtleSkgPT4gT2JqZWN0LmFzc2lnbihyZXN1bHQsIHsgW2tleV06IG9ba2V5XSB9KSxcbiAgICB7fVxuICApO1xufTtcblxuLyoqXG4gKiBFeGNsdWRlIGtleXMgZnJvbSBhbiBvYmplY3RcbiAqIEBwYXJhbSB7b2JqZWN0fSBvXG4gKiBAcGFyYW0gey4uLnN0cmluZ30ga2V5cyBUaGUga2V5cyB0byBleGNsdWRlIGZyb20gdGhlIG9iamVjdFxuICogQHJldHVybiB7b2JqZWN0fSBBbiBvYmplY3QgY29udGFpbmluZyBhbGwga2V5cyBleGNlcHQgZXhjbHVkZWQga2V5c1xuICovXG5jb25zdCBleGNsdWRlID0gKG8sIC4uLmtleXMpID0+IHtcbiAgcmV0dXJuIE9iamVjdC5mcm9tRW50cmllcyhcbiAgICBPYmplY3QuZW50cmllcyhvKS5maWx0ZXIoKFtrZXldKSA9PiAha2V5cy5pbmNsdWRlcyhrZXkpKVxuICApO1xufTtcblxuaWYgKHR5cGVvZiBtb2R1bGUgIT09ICd1bmRlZmluZWQnKSB7XG4gIG1vZHVsZS5leHBvcnRzID0ge1xuICAgIGZsb2F0RXF1YWxzLFxuICAgIGNsYW1wLFxuICAgIGZyYWMsXG4gICAgcm91bmQsXG4gICAgbGVycCxcbiAgICB1bmxlcnAsXG4gICAgYmxlcnAsXG4gICAgcmVtYXAsXG4gICAgc21vb3Roc3RlcCxcbiAgICByYWRpYW5zLFxuICAgIGRlZ3JlZXMsXG4gICAgcmFuZG9tQmV0d2VlbixcbiAgICByYW5kb21JbnRCZXR3ZWVuLFxuICAgIGNsdFJhbmRvbSxcbiAgICBjbHRSYW5kb21JbnQsXG4gICAgd2VpZ2h0ZWRSYW5kb20sXG4gICAgbGVycEFycmF5LFxuICAgIGRvdCxcbiAgICBmYWN0b3JpYWwsXG4gICAgbnByLFxuICAgIG5jcixcbiAgICBjb21iaW5hdGlvbnMsXG4gICAgY2FydGVzaWFuLFxuICAgIHRpbWVzLFxuICAgIHJhbmdlLFxuICAgIHppcCxcbiAgICBhdCxcbiAgICBwZWVrLFxuICAgIGNodW5rLFxuICAgIHNodWZmbGUsXG4gICAgZmxhdCxcbiAgICB1bmZsYXQsXG4gICAgc3BsaXQsXG4gICAgcGx1Y2ssXG4gICAgZXhjbHVkZSxcbiAgfTtcbn1cbiIsImNvbnN0IHsgdGltZXMsIGNodW5rLCBkb3QgfSA9IHJlcXVpcmUoJ0BiYXNlbWVudHVuaXZlcnNlL3V0aWxzJyk7XG5cbi8qKlxuICogQG92ZXJ2aWV3IEEgc21hbGwgdmVjdG9yIGFuZCBtYXRyaXggbGlicmFyeVxuICogQGF1dGhvciBHb3Jkb24gTGFycmlnYW5cbiAqL1xuXG4vKipcbiAqIEEgMmQgdmVjdG9yXG4gKiBAdHlwZWRlZiB7T2JqZWN0fSB2ZWNcbiAqIEBwcm9wZXJ0eSB7bnVtYmVyfSB4IFRoZSB4IGNvbXBvbmVudCBvZiB0aGUgdmVjdG9yXG4gKiBAcHJvcGVydHkge251bWJlcn0geSBUaGUgeSBjb21wb25lbnQgb2YgdGhlIHZlY3RvclxuICovXG5cbi8qKlxuICogQ3JlYXRlIGEgbmV3IHZlY3RvclxuICogQHBhcmFtIHtudW1iZXJ8dmVjfSBbeF0gVGhlIHggY29tcG9uZW50IG9mIHRoZSB2ZWN0b3IsIG9yIGEgdmVjdG9yIHRvIGNvcHlcbiAqIEBwYXJhbSB7bnVtYmVyfSBbeV0gVGhlIHkgY29tcG9uZW50IG9mIHRoZSB2ZWN0b3JcbiAqIEByZXR1cm4ge3ZlY30gQSBuZXcgdmVjdG9yXG4gKiBAZXhhbXBsZSA8Y2FwdGlvbj5WYXJpb3VzIHdheXMgdG8gaW5pdGlhbGlzZSBhIHZlY3RvcjwvY2FwdGlvbj5cbiAqIGxldCBhID0gdmVjKDMsIDIpOyAgLy8gKDMsIDIpXG4gKiBsZXQgYiA9IHZlYyg0KTsgICAgIC8vICg0LCA0KVxuICogbGV0IGMgPSB2ZWMoYSk7ICAgICAvLyAoMywgMilcbiAqIGxldCBkID0gdmVjKCk7ICAgICAgLy8gKDAsIDApXG4gKi9cbmNvbnN0IHZlYyA9ICh4LCB5KSA9PiAoIXggJiYgIXkgP1xuICB7IHg6IDAsIHk6IDAgfSA6ICh0eXBlb2YgeCA9PT0gJ29iamVjdCcgP1xuICAgIHsgeDogeC54IHx8IDAsIHk6IHgueSB8fCAwIH0gOiAoeSA9PT0gbnVsbCB8fCB5ID09PSB1bmRlZmluZWQgP1xuICAgICAgeyB4OiB4LCB5OiB4IH0gOiB7IHg6IHgsIHk6IHkgfSlcbiAgKVxuKTtcblxuLyoqXG4gKiBHZXQgdGhlIGNvbXBvbmVudHMgb2YgYSB2ZWN0b3IgYXMgYW4gYXJyYXlcbiAqIEBwYXJhbSB7dmVjfSBhIFRoZSB2ZWN0b3IgdG8gZ2V0IGNvbXBvbmVudHMgZnJvbVxuICogQHJldHVybiB7QXJyYXk8bnVtYmVyPn0gVGhlIHZlY3RvciBjb21wb25lbnRzIGFzIGFuIGFycmF5XG4gKi9cbnZlYy5jb21wb25lbnRzID0gYSA9PiBbYS54LCBhLnldO1xuXG4vKipcbiAqIFJldHVybiBhIHVuaXQgdmVjdG9yICgxLCAwKVxuICogQHJldHVybiB7dmVjfSBBIHVuaXQgdmVjdG9yICgxLCAwKVxuICovXG52ZWMudXggPSAoKSA9PiB2ZWMoMSwgMCk7XG5cbi8qKlxuICogUmV0dXJuIGEgdW5pdCB2ZWN0b3IgKDAsIDEpXG4gKiBAcmV0dXJuIHt2ZWN9IEEgdW5pdCB2ZWN0b3IgKDAsIDEpXG4gKi9cbnZlYy51eSA9ICgpID0+IHZlYygwLCAxKTtcblxuLyoqXG4gKiBBZGQgdmVjdG9yc1xuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEBwYXJhbSB7dmVjfSBiIFZlY3RvciBiXG4gKiBAcmV0dXJuIHt2ZWN9IGEgKyBiXG4gKi9cbnZlYy5hZGQgPSAoYSwgYikgPT4gKHsgeDogYS54ICsgYi54LCB5OiBhLnkgKyBiLnkgfSk7XG5cbi8qKlxuICogU2NhbGUgYSB2ZWN0b3JcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcGFyYW0ge251bWJlcn0gYiBTY2FsYXIgYlxuICogQHJldHVybiB7dmVjfSBhICogYlxuICovXG52ZWMubXVsID0gKGEsIGIpID0+ICh7IHg6IGEueCAqIGIsIHk6IGEueSAqIGIgfSk7XG5cbi8qKlxuICogU3VidHJhY3QgdmVjdG9yc1xuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEBwYXJhbSB7dmVjfSBiIFZlY3RvciBiXG4gKiBAcmV0dXJuIHt2ZWN9IGEgLSBiXG4gKi9cbnZlYy5zdWIgPSAoYSwgYikgPT4gKHsgeDogYS54IC0gYi54LCB5OiBhLnkgLSBiLnkgfSk7XG5cbi8qKlxuICogR2V0IHRoZSBsZW5ndGggb2YgYSB2ZWN0b3JcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcmV0dXJuIHtudW1iZXJ9IHxhfFxuICovXG52ZWMubGVuID0gYSA9PiBNYXRoLnNxcnQoYS54ICogYS54ICsgYS55ICogYS55KTtcblxuLyoqXG4gKiBHZXQgdGhlIGxlbmd0aCBvZiBhIHZlY3RvciB1c2luZyB0YXhpY2FiIGdlb21ldHJ5XG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHJldHVybiB7bnVtYmVyfSB8YXxcbiAqL1xudmVjLm1hbmhhdHRhbiA9IGEgPT4gTWF0aC5hYnMoYS54KSArIE1hdGguYWJzKGEueSk7XG5cbi8qKlxuICogTm9ybWFsaXNlIGEgdmVjdG9yXG4gKiBAcGFyYW0ge3ZlY30gYSBUaGUgdmVjdG9yIHRvIG5vcm1hbGlzZVxuICogQHJldHVybiB7dmVjfSBeYVxuICovXG52ZWMubm9yID0gYSA9PiB7XG4gIGxldCBsZW4gPSB2ZWMubGVuKGEpO1xuICByZXR1cm4gbGVuID8geyB4OiBhLnggLyBsZW4sIHk6IGEueSAvIGxlbiB9IDogdmVjKCk7XG59O1xuXG4vKipcbiAqIEdldCBhIGRvdCBwcm9kdWN0IG9mIHZlY3RvcnNcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcGFyYW0ge3ZlY30gYiBWZWN0b3IgYlxuICogQHJldHVybiB7bnVtYmVyfSBhIOKImSBiXG4gKi9cbnZlYy5kb3QgPSAoYSwgYikgPT4gYS54ICogYi54ICsgYS55ICogYi55O1xuXG4vKipcbiAqIFJvdGF0ZSBhIHZlY3RvciBieSByIHJhZGlhbnNcbiAqIEBwYXJhbSB7dmVjfSBhIFRoZSB2ZWN0b3IgdG8gcm90YXRlXG4gKiBAcGFyYW0ge251bWJlcn0gciBUaGUgYW5nbGUgdG8gcm90YXRlIGJ5LCBtZWFzdXJlZCBpbiByYWRpYW5zXG4gKiBAcmV0dXJuIHt2ZWN9IEEgcm90YXRlZCB2ZWN0b3JcbiAqL1xudmVjLnJvdCA9IChhLCByKSA9PiB7XG4gIGxldCBzID0gTWF0aC5zaW4ociksXG4gICAgYyA9IE1hdGguY29zKHIpO1xuICByZXR1cm4geyB4OiBjICogYS54IC0gcyAqIGEueSwgeTogcyAqIGEueCArIGMgKiBhLnkgfTtcbn1cblxuLyoqXG4gKiBDaGVjayBpZiB0d28gdmVjdG9ycyBhcmUgZXF1YWxcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcGFyYW0ge3ZlY30gYiBWZWN0b3IgYlxuICogQHJldHVybiB7Ym9vbGVhbn0gVHJ1ZSBpZiB2ZWN0b3JzIGEgYW5kIGIgYXJlIGVxdWFsLCBmYWxzZSBvdGhlcndpc2VcbiAqL1xudmVjLmVxID0gKGEsIGIpID0+IGEueCA9PT0gYi54ICYmIGEueSA9PT0gYi55O1xuXG4vKipcbiAqIEdldCB0aGUgYW5nbGUgb2YgYSB2ZWN0b3JcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcmV0dXJuIHtudW1iZXJ9IFRoZSBhbmdsZSBvZiB2ZWN0b3IgYSBpbiByYWRpYW5zXG4gKi9cbnZlYy5yYWQgPSBhID0+IE1hdGguYXRhbjIoYS55LCBhLngpO1xuXG4vKipcbiAqIENvcHkgYSB2ZWN0b3JcbiAqIEBwYXJhbSB7dmVjfSBhIFRoZSB2ZWN0b3IgdG8gY29weVxuICogQHJldHVybiB7dmVjfSBBIGNvcHkgb2YgdmVjdG9yIGFcbiAqL1xudmVjLmNweSA9IGEgPT4gdmVjKGEpO1xuXG4vKipcbiAqIEEgZnVuY3Rpb24gdG8gY2FsbCBvbiBlYWNoIGNvbXBvbmVudCBvZiBhIHZlY3RvclxuICogQGNhbGxiYWNrIHZlY3Rvck1hcENhbGxiYWNrXG4gKiBAcGFyYW0ge251bWJlcn0gdmFsdWUgVGhlIGNvbXBvbmVudCB2YWx1ZVxuICogQHBhcmFtIHsneCcgfCAneSd9IGxhYmVsIFRoZSBjb21wb25lbnQgbGFiZWwgKHggb3IgeSlcbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIG1hcHBlZCBjb21wb25lbnRcbiAqL1xuXG4vKipcbiAqIENhbGwgYSBmdW5jdGlvbiBvbiBlYWNoIGNvbXBvbmVudCBvZiBhIHZlY3RvciBhbmQgYnVpbGQgYSBuZXcgdmVjdG9yIGZyb20gdGhlIHJlc3VsdHNcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcGFyYW0ge3ZlY3Rvck1hcENhbGxiYWNrfSBmIFRoZSBmdW5jdGlvbiB0byBjYWxsIG9uIGVhY2ggY29tcG9uZW50IG9mIHRoZSB2ZWN0b3JcbiAqIEByZXR1cm4ge3ZlY30gVmVjdG9yIGEgbWFwcGVkIHRocm91Z2ggZlxuICovXG52ZWMubWFwID0gKGEsIGYpID0+ICh7IHg6IGYoYS54LCAneCcpLCB5OiBmKGEueSwgJ3knKSB9KTtcblxuLyoqXG4gKiBDb252ZXJ0IGEgdmVjdG9yIGludG8gYSBzdHJpbmdcbiAqIEBwYXJhbSB7dmVjfSBhIFRoZSB2ZWN0b3IgdG8gY29udmVydFxuICogQHBhcmFtIHtzdHJpbmd9IFtzPScsICddIFRoZSBzZXBhcmF0b3Igc3RyaW5nXG4gKiBAcmV0dXJuIHtzdHJpbmd9IEEgc3RyaW5nIHJlcHJlc2VudGF0aW9uIG9mIHRoZSB2ZWN0b3JcbiAqL1xudmVjLnN0ciA9IChhLCBzID0gJywgJykgPT4gYCR7YS54fSR7c30ke2EueX1gO1xuXG4vKipcbiAqIEEgbWF0cml4XG4gKiBAdHlwZWRlZiB7T2JqZWN0fSBtYXRcbiAqIEBwcm9wZXJ0eSB7bnVtYmVyfSBtIFRoZSBudW1iZXIgb2Ygcm93cyBpbiB0aGUgbWF0cml4XG4gKiBAcHJvcGVydHkge251bWJlcn0gbiBUaGUgbnVtYmVyIG9mIGNvbHVtbnMgaW4gdGhlIG1hdHJpeFxuICogQHByb3BlcnR5IHtBcnJheTxudW1iZXI+fSBlbnRyaWVzIFRoZSBtYXRyaXggdmFsdWVzXG4gKi9cblxuLyoqXG4gKiBDcmVhdGUgYSBuZXcgbWF0cml4XG4gKiBAcGFyYW0ge251bWJlcn0gW209NF0gVGhlIG51bWJlciBvZiByb3dzXG4gKiBAcGFyYW0ge251bWJlcn0gW249NF0gVGhlIG51bWJlciBvZiBjb2x1bW5zXG4gKiBAcGFyYW0ge0FycmF5PG51bWJlcj59IFtlbnRyaWVzPVtdXSBNYXRyaXggdmFsdWVzIGluIHJlYWRpbmcgb3JkZXJcbiAqIEByZXR1cm4ge21hdH0gQSBuZXcgbWF0cml4XG4gKi9cbmNvbnN0IG1hdCA9IChtID0gNCwgbiA9IDQsIGVudHJpZXMgPSBbXSkgPT4gKHtcbiAgbSwgbixcbiAgZW50cmllczogZW50cmllcy5jb25jYXQoQXJyYXkobSAqIG4pLmZpbGwoMCkpLnNsaWNlKDAsIG0gKiBuKVxufSk7XG5cbi8qKlxuICogR2V0IGFuIGlkZW50aXR5IG1hdHJpeCBvZiBzaXplIG5cbiAqIEBwYXJhbSB7bnVtYmVyfSBuIFRoZSBzaXplIG9mIHRoZSBtYXRyaXhcbiAqIEByZXR1cm4ge21hdH0gQW4gaWRlbnRpdHkgbWF0cml4XG4gKi9cbm1hdC5pZGVudGl0eSA9IG4gPT4gbWF0KG4sIG4sIEFycmF5KG4gKiBuKS5maWxsKDApLm1hcCgodiwgaSkgPT4gKyhNYXRoLmZsb29yKGkgLyBuKSA9PT0gaSAlIG4pKSk7XG5cbi8qKlxuICogR2V0IGFuIGVudHJ5IGZyb20gYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgcm93IG9mZnNldFxuICogQHBhcmFtIHtudW1iZXJ9IGogVGhlIGNvbHVtbiBvZmZzZXRcbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIHZhbHVlIGF0IHBvc2l0aW9uIChpLCBqKSBpbiBtYXRyaXggYVxuICovXG5tYXQuZ2V0ID0gKGEsIGksIGopID0+IGEuZW50cmllc1soaiAtIDEpICsgKGkgLSAxKSAqIGEubl07XG5cbi8qKlxuICogU2V0IGFuIGVudHJ5IG9mIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIHJvdyBvZmZzZXRcbiAqIEBwYXJhbSB7bnVtYmVyfSBqIFRoZSBjb2x1bW4gb2Zmc2V0XG4gKiBAcGFyYW0ge251bWJlcn0gdiBUaGUgdmFsdWUgdG8gc2V0IGluIG1hdHJpeCBhXG4gKi9cbm1hdC5zZXQgPSAoYSwgaSwgaiwgdikgPT4geyBhLmVudHJpZXNbKGogLSAxKSArIChpIC0gMSkgKiBhLm5dID0gdjsgfTtcblxuLyoqXG4gKiBHZXQgYSByb3cgZnJvbSBhIG1hdHJpeCBhcyBhbiBhcnJheVxuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBtIFRoZSByb3cgb2Zmc2V0XG4gKiBAcmV0dXJuIHtBcnJheTxudW1iZXI+fSBSb3cgbSBmcm9tIG1hdHJpeCBhXG4gKi9cbm1hdC5yb3cgPSAoYSwgbSkgPT4ge1xuICBjb25zdCBzID0gKG0gLSAxKSAqIGEubjtcbiAgcmV0dXJuIGEuZW50cmllcy5zbGljZShzLCBzICsgYS5uKTtcbn07XG5cbi8qKlxuICogR2V0IGEgY29sdW1uIGZyb20gYSBtYXRyaXggYXMgYW4gYXJyYXlcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge251bWJlcn0gbiBUaGUgY29sdW1uIG9mZnNldFxuICogQHJldHVybiB7QXJyYXk8bnVtYmVyPn0gQ29sdW1uIG4gZnJvbSBtYXRyaXggYVxuICovXG5tYXQuY29sID0gKGEsIG4pID0+IHRpbWVzKGkgPT4gbWF0LmdldChhLCAoaSArIDEpLCBuKSwgYS5tKTtcblxuLyoqXG4gKiBBZGQgbWF0cmljZXNcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge21hdH0gYiBNYXRyaXggYlxuICogQHJldHVybiB7bWF0fSBhICsgYlxuICovXG5tYXQuYWRkID0gKGEsIGIpID0+IGEubSA9PT0gYi5tICYmIGEubiA9PT0gYi5uICYmIG1hdC5tYXAoYSwgKHYsIGkpID0+IHYgKyBiLmVudHJpZXNbaV0pO1xuXG4vKipcbiAqIFN1YnRyYWN0IG1hdHJpY2VzXG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHttYXR9IGIgTWF0cml4IGJcbiAqIEByZXR1cm4ge21hdH0gYSAtIGJcbiAqL1xubWF0LnN1YiA9IChhLCBiKSA9PiBhLm0gPT09IGIubSAmJiBhLm4gPT09IGIubiAmJiBtYXQubWFwKGEsICh2LCBpKSA9PiB2IC0gYi5lbnRyaWVzW2ldKTtcblxuLyoqXG4gKiBNdWx0aXBseSBtYXRyaWNlc1xuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bWF0fSBiIE1hdHJpeCBiXG4gKiBAcmV0dXJuIHttYXR8Ym9vbGVhbn0gYWIgb3IgZmFsc2UgaWYgdGhlIG1hdHJpY2VzIGNhbm5vdCBiZSBtdWx0aXBsaWVkXG4gKi9cbm1hdC5tdWwgPSAoYSwgYikgPT4ge1xuICBpZiAoYS5uICE9PSBiLm0pIHsgcmV0dXJuIGZhbHNlOyB9XG4gIGNvbnN0IHJlc3VsdCA9IG1hdChhLm0sIGIubik7XG4gIGZvciAobGV0IGkgPSAxOyBpIDw9IGEubTsgaSsrKSB7XG4gICAgZm9yIChsZXQgaiA9IDE7IGogPD0gYi5uOyBqKyspIHtcbiAgICAgIG1hdC5zZXQocmVzdWx0LCBpLCBqLCBkb3QobWF0LnJvdyhhLCBpKSwgbWF0LmNvbChiLCBqKSkpO1xuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzdWx0O1xufTtcblxuLyoqXG4gKiBTY2FsZSBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBiIFNjYWxhciBiXG4gKiBAcmV0dXJuIHttYXR9IGEgKiBiXG4gKi9cbm1hdC5zY2FsZSA9IChhLCBiKSA9PiBtYXQubWFwKGEsIHYgPT4gdiAqIGIpO1xuXG4vKipcbiAqIFRyYW5zcG9zZSBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgVGhlIG1hdHJpeCB0byB0cmFuc3Bvc2VcbiAqIEByZXR1cm4ge21hdH0gQSB0cmFuc3Bvc2VkIG1hdHJpeFxuICovXG5tYXQudHJhbnMgPSBhID0+IG1hdChhLm4sIGEubSwgdGltZXMoaSA9PiBtYXQuY29sKGEsIChpICsgMSkpLCBhLm4pLmZsYXQoKSk7XG5cbi8qKlxuICogR2V0IHRoZSBtaW5vciBvZiBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSByb3cgb2Zmc2V0XG4gKiBAcGFyYW0ge251bWJlcn0gaiBUaGUgY29sdW1uIG9mZnNldFxuICogQHJldHVybiB7bWF0fGJvb2xlYW59IFRoZSAoaSwgaikgbWlub3Igb2YgbWF0cml4IGEgb3IgZmFsc2UgaWYgdGhlIG1hdHJpeCBpcyBub3Qgc3F1YXJlXG4gKi9cbm1hdC5taW5vciA9IChhLCBpLCBqKSA9PiB7XG4gIGlmIChhLm0gIT09IGEubikgeyByZXR1cm4gZmFsc2U7IH1cbiAgY29uc3QgZW50cmllcyA9IFtdO1xuICBmb3IgKGxldCBpaSA9IDE7IGlpIDw9IGEubTsgaWkrKykge1xuICAgIGlmIChpaSA9PT0gaSkgeyBjb250aW51ZTsgfVxuICAgIGZvciAobGV0IGpqID0gMTsgamogPD0gYS5uOyBqaisrKSB7XG4gICAgICBpZiAoamogPT09IGopIHsgY29udGludWU7IH1cbiAgICAgIGVudHJpZXMucHVzaChtYXQuZ2V0KGEsIGlpLCBqaikpO1xuICAgIH1cbiAgfVxuICByZXR1cm4gbWF0KGEubSAtIDEsIGEubiAtIDEsIGVudHJpZXMpO1xufTtcblxuLyoqXG4gKiBHZXQgdGhlIGRldGVybWluYW50IG9mIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHJldHVybiB7bnVtYmVyfGJvb2xlYW59IHxhfCBvciBmYWxzZSBpZiB0aGUgbWF0cml4IGlzIG5vdCBzcXVhcmVcbiAqL1xubWF0LmRldCA9IGEgPT4ge1xuICBpZiAoYS5tICE9PSBhLm4pIHsgcmV0dXJuIGZhbHNlOyB9XG4gIGlmIChhLm0gPT09IDEpIHtcbiAgICByZXR1cm4gYS5lbnRyaWVzWzBdO1xuICB9XG4gIGlmIChhLm0gPT09IDIpIHtcbiAgICByZXR1cm4gYS5lbnRyaWVzWzBdICogYS5lbnRyaWVzWzNdIC0gYS5lbnRyaWVzWzFdICogYS5lbnRyaWVzWzJdO1xuICB9XG4gIGxldCB0b3RhbCA9IDAsIHNpZ24gPSAxO1xuICBmb3IgKGxldCBqID0gMTsgaiA8PSBhLm47IGorKykge1xuICAgIHRvdGFsICs9IHNpZ24gKiBhLmVudHJpZXNbaiAtIDFdICogbWF0LmRldChtYXQubWlub3IoYSwgMSwgaikpO1xuICAgIHNpZ24gKj0gLTE7XG4gIH1cbiAgcmV0dXJuIHRvdGFsO1xufTtcblxuLyoqXG4gKiBOb3JtYWxpc2UgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIFRoZSBtYXRyaXggdG8gbm9ybWFsaXNlXG4gKiBAcmV0dXJuIHttYXR8Ym9vbGVhbn0gXmEgb3IgZmFsc2UgaWYgdGhlIG1hdHJpeCBpcyBub3Qgc3F1YXJlXG4gKi9cbm1hdC5ub3IgPSBhID0+IHtcbiAgaWYgKGEubSAhPT0gYS5uKSB7IHJldHVybiBmYWxzZTsgfVxuICBjb25zdCBkID0gbWF0LmRldChhKTtcbiAgcmV0dXJuIG1hdC5tYXAoYSwgaSA9PiBpICogZCk7XG59O1xuXG4vKipcbiAqIEdldCB0aGUgYWRqdWdhdGUgb2YgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIFRoZSBtYXRyaXggZnJvbSB3aGljaCB0byBnZXQgdGhlIGFkanVnYXRlXG4gKiBAcmV0dXJuIHttYXR9IFRoZSBhZGp1Z2F0ZSBvZiBhXG4gKi9cbm1hdC5hZGogPSBhID0+IHtcbiAgY29uc3QgbWlub3JzID0gbWF0KGEubSwgYS5uKTtcbiAgZm9yIChsZXQgaSA9IDE7IGkgPD0gYS5tOyBpKyspIHtcbiAgICBmb3IgKGxldCBqID0gMTsgaiA8PSBhLm47IGorKykge1xuICAgICAgbWF0LnNldChtaW5vcnMsIGksIGosIG1hdC5kZXQobWF0Lm1pbm9yKGEsIGksIGopKSk7XG4gICAgfVxuICB9XG4gIGNvbnN0IGNvZmFjdG9ycyA9IG1hdC5tYXAobWlub3JzLCAodiwgaSkgPT4gdiAqIChpICUgMiA/IC0xIDogMSkpO1xuICByZXR1cm4gbWF0LnRyYW5zKGNvZmFjdG9ycyk7XG59O1xuXG4vKipcbiAqIEdldCB0aGUgaW52ZXJzZSBvZiBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgVGhlIG1hdHJpeCB0byBpbnZlcnRcbiAqIEByZXR1cm4ge21hdHxib29sZWFufSBhXi0xIG9yIGZhbHNlIGlmIHRoZSBtYXRyaXggaGFzIG5vIGludmVyc2VcbiAqL1xubWF0LmludiA9IGEgPT4ge1xuICBpZiAoYS5tICE9PSBhLm4pIHsgcmV0dXJuIGZhbHNlOyB9XG4gIGNvbnN0IGQgPSBtYXQuZGV0KGEpO1xuICBpZiAoZCA9PT0gMCkgeyByZXR1cm4gZmFsc2U7IH1cbiAgcmV0dXJuIG1hdC5zY2FsZShtYXQuYWRqKGEpLCAxIC8gZCk7XG59O1xuXG4vKipcbiAqIENoZWNrIGlmIHR3byBtYXRyaWNlcyBhcmUgZXF1YWxcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge21hdH0gYiBNYXRyaXggYlxuICogQHJldHVybiB7Ym9vbGVhbn0gVHJ1ZSBpZiBtYXRyaWNlcyBhIGFuZCBiIGFyZSBpZGVudGljYWwsIGZhbHNlIG90aGVyd2lzZVxuICovXG5tYXQuZXEgPSAoYSwgYikgPT4gYS5tID09PSBiLm0gJiYgYS5uID09PSBiLm4gJiYgbWF0LnN0cihhKSA9PT0gbWF0LnN0cihiKTtcblxuLyoqXG4gKiBDb3B5IGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBUaGUgbWF0cml4IHRvIGNvcHlcbiAqIEByZXR1cm4ge21hdH0gQSBjb3B5IG9mIG1hdHJpeCBhXG4gKi9cbm1hdC5jcHkgPSBhID0+IG1hdChhLm0sIGEubiwgWy4uLmEuZW50cmllc10pO1xuXG4vKipcbiAqIEEgZnVuY3Rpb24gdG8gY2FsbCBvbiBlYWNoIGVudHJ5IG9mIGEgbWF0cml4XG4gKiBAY2FsbGJhY2sgbWF0cml4TWFwQ2FsbGJhY2tcbiAqIEBwYXJhbSB7bnVtYmVyfSB2YWx1ZSBUaGUgZW50cnkgdmFsdWVcbiAqIEBwYXJhbSB7bnVtYmVyfSBpbmRleCBUaGUgZW50cnkgaW5kZXhcbiAqIEBwYXJhbSB7QXJyYXk8bnVtYmVyPn0gZW50cmllcyBUaGUgYXJyYXkgb2YgbWF0cml4IGVudHJpZXNcbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIG1hcHBlZCBlbnRyeVxuICovXG5cbi8qKlxuICogQ2FsbCBhIGZ1bmN0aW9uIG9uIGVhY2ggZW50cnkgb2YgYSBtYXRyaXggYW5kIGJ1aWxkIGEgbmV3IG1hdHJpeCBmcm9tIHRoZSByZXN1bHRzXG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHttYXRyaXhNYXBDYWxsYmFja30gZiBUaGUgZnVuY3Rpb24gdG8gY2FsbCBvbiBlYWNoIGVudHJ5IG9mIHRoZSBtYXRyaXhcbiAqIEByZXR1cm4ge21hdH0gTWF0cml4IGEgbWFwcGVkIHRocm91Z2ggZlxuICovXG5tYXQubWFwID0gKGEsIGYpID0+IG1hdChhLm0sIGEubiwgYS5lbnRyaWVzLm1hcChmKSk7XG5cbi8qKlxuICogQ29udmVydCBhIG1hdHJpeCBpbnRvIGEgc3RyaW5nXG4gKiBAcGFyYW0ge21hdH0gYSBUaGUgbWF0cml4IHRvIGNvbnZlcnRcbiAqIEBwYXJhbSB7c3RyaW5nfSBbbXM9JywgJ10gVGhlIHNlcGFyYXRvciBzdHJpbmcgZm9yIGNvbHVtbnNcbiAqIEBwYXJhbSB7c3RyaW5nfSBbbnM9J1xcbiddIFRoZSBzZXBhcmF0b3Igc3RyaW5nIGZvciByb3dzXG4gKiBAcmV0dXJuIHtzdHJpbmd9IEEgc3RyaW5nIHJlcHJlc2VudGF0aW9uIG9mIHRoZSBtYXRyaXhcbiAqL1xubWF0LnN0ciA9IChhLCBtcyA9ICcsICcsIG5zID0gJ1xcbicpID0+IGNodW5rKGEuZW50cmllcywgYS5uKS5tYXAociA9PiByLmpvaW4obXMpKS5qb2luKG5zKTtcblxuaWYgKHR5cGVvZiBtb2R1bGUgIT09ICd1bmRlZmluZWQnKSB7XG4gIG1vZHVsZS5leHBvcnRzID0geyB2ZWMsIG1hdCB9O1xufVxuIiwiIWZ1bmN0aW9uKGcsYyl7dHlwZW9mIGV4cG9ydHM9PVwib2JqZWN0XCImJnR5cGVvZiBtb2R1bGUhPVwidW5kZWZpbmVkXCI/YyhleHBvcnRzKTp0eXBlb2YgZGVmaW5lPT1cImZ1bmN0aW9uXCImJmRlZmluZS5hbWQ/ZGVmaW5lKFtcImV4cG9ydHNcIl0sYyk6YygoZz1nfHxzZWxmKS5scnVfbWFwPWcubHJ1X21hcHx8e30pfSh0aGlzLGZ1bmN0aW9uKGcpe2NvbnN0IGM9U3ltYm9sKFwibmV3ZXJcIiksZT1TeW1ib2woXCJvbGRlclwiKTtjbGFzcyBue2NvbnN0cnVjdG9yKGEsYil7dHlwZW9mIGEhPT1cIm51bWJlclwiJiYoYj1hLGE9MCksdGhpcy5zaXplPTAsdGhpcy5saW1pdD1hLHRoaXMub2xkZXN0PXRoaXMubmV3ZXN0PXZvaWQgMCx0aGlzLl9rZXltYXA9bmV3IE1hcCgpLGImJih0aGlzLmFzc2lnbihiKSxhPDEmJih0aGlzLmxpbWl0PXRoaXMuc2l6ZSkpfV9tYXJrRW50cnlBc1VzZWQoYSl7aWYoYT09PXRoaXMubmV3ZXN0KXJldHVybjthW2NdJiYoYT09PXRoaXMub2xkZXN0JiYodGhpcy5vbGRlc3Q9YVtjXSksYVtjXVtlXT1hW2VdKSxhW2VdJiYoYVtlXVtjXT1hW2NdKSxhW2NdPXZvaWQgMCxhW2VdPXRoaXMubmV3ZXN0LHRoaXMubmV3ZXN0JiYodGhpcy5uZXdlc3RbY109YSksdGhpcy5uZXdlc3Q9YX1hc3NpZ24oYSl7bGV0IGIsZD10aGlzLmxpbWl0fHxOdW1iZXIuTUFYX1ZBTFVFO3RoaXMuX2tleW1hcC5jbGVhcigpO2xldCBtPWFbU3ltYm9sLml0ZXJhdG9yXSgpO2ZvcihsZXQgaD1tLm5leHQoKTshaC5kb25lO2g9bS5uZXh0KCkpe2xldCBmPW5ldyBsKGgudmFsdWVbMF0saC52YWx1ZVsxXSk7dGhpcy5fa2V5bWFwLnNldChmLmtleSxmKSxiPyhiW2NdPWYsZltlXT1iKTp0aGlzLm9sZGVzdD1mLGI9ZjtpZihkLS09PTApdGhyb3cgbmV3IEVycm9yKFwib3ZlcmZsb3dcIil9dGhpcy5uZXdlc3Q9Yix0aGlzLnNpemU9dGhpcy5fa2V5bWFwLnNpemV9Z2V0KGEpe3ZhciBiPXRoaXMuX2tleW1hcC5nZXQoYSk7cmV0dXJuIGI/KHRoaXMuX21hcmtFbnRyeUFzVXNlZChiKSxiLnZhbHVlKTp2b2lkIDB9c2V0KGEsYil7dmFyIGQ9dGhpcy5fa2V5bWFwLmdldChhKTtyZXR1cm4gZD8oZC52YWx1ZT1iLHRoaXMuX21hcmtFbnRyeUFzVXNlZChkKSx0aGlzKToodGhpcy5fa2V5bWFwLnNldChhLGQ9bmV3IGwoYSxiKSksdGhpcy5uZXdlc3Q/KHRoaXMubmV3ZXN0W2NdPWQsZFtlXT10aGlzLm5ld2VzdCk6dGhpcy5vbGRlc3Q9ZCx0aGlzLm5ld2VzdD1kLCsrdGhpcy5zaXplLHRoaXMuc2l6ZT50aGlzLmxpbWl0JiZ0aGlzLnNoaWZ0KCksdGhpcyl9c2hpZnQoKXt2YXIgYT10aGlzLm9sZGVzdDtpZihhKXJldHVybiB0aGlzLm9sZGVzdFtjXT8odGhpcy5vbGRlc3Q9dGhpcy5vbGRlc3RbY10sdGhpcy5vbGRlc3RbZV09dm9pZCAwKToodGhpcy5vbGRlc3Q9dm9pZCAwLHRoaXMubmV3ZXN0PXZvaWQgMCksYVtjXT1hW2VdPXZvaWQgMCx0aGlzLl9rZXltYXAuZGVsZXRlKGEua2V5KSwtLXRoaXMuc2l6ZSxbYS5rZXksYS52YWx1ZV19ZmluZChhKXtsZXQgYj10aGlzLl9rZXltYXAuZ2V0KGEpO3JldHVybiBiP2IudmFsdWU6dm9pZCAwfWhhcyhhKXtyZXR1cm4gdGhpcy5fa2V5bWFwLmhhcyhhKX1kZWxldGUoYSl7dmFyIGI9dGhpcy5fa2V5bWFwLmdldChhKTtyZXR1cm4gYj8odGhpcy5fa2V5bWFwLmRlbGV0ZShiLmtleSksYltjXSYmYltlXT8oYltlXVtjXT1iW2NdLGJbY11bZV09YltlXSk6YltjXT8oYltjXVtlXT12b2lkIDAsdGhpcy5vbGRlc3Q9YltjXSk6YltlXT8oYltlXVtjXT12b2lkIDAsdGhpcy5uZXdlc3Q9YltlXSk6dGhpcy5vbGRlc3Q9dGhpcy5uZXdlc3Q9dm9pZCAwLHRoaXMuc2l6ZS0tLGIudmFsdWUpOnZvaWQgMH1jbGVhcigpe3RoaXMub2xkZXN0PXRoaXMubmV3ZXN0PXZvaWQgMCx0aGlzLnNpemU9MCx0aGlzLl9rZXltYXAuY2xlYXIoKX1rZXlzKCl7cmV0dXJuIG5ldyBqKHRoaXMub2xkZXN0KX12YWx1ZXMoKXtyZXR1cm4gbmV3IGsodGhpcy5vbGRlc3QpfWVudHJpZXMoKXtyZXR1cm4gdGhpc31bU3ltYm9sLml0ZXJhdG9yXSgpe3JldHVybiBuZXcgaSh0aGlzLm9sZGVzdCl9Zm9yRWFjaChhLGIpe3R5cGVvZiBiIT09XCJvYmplY3RcIiYmKGI9dGhpcyk7bGV0IGQ9dGhpcy5vbGRlc3Q7Zm9yKDtkOylhLmNhbGwoYixkLnZhbHVlLGQua2V5LHRoaXMpLGQ9ZFtjXX10b0pTT04oKXtmb3IodmFyIGE9bmV3IEFycmF5KHRoaXMuc2l6ZSksYj0wLGQ9dGhpcy5vbGRlc3Q7ZDspYVtiKytdPXtrZXk6ZC5rZXksdmFsdWU6ZC52YWx1ZX0sZD1kW2NdO3JldHVybiBhfXRvU3RyaW5nKCl7Zm9yKHZhciBhPVwiXCIsYj10aGlzLm9sZGVzdDtiOylhKz1TdHJpbmcoYi5rZXkpK1wiOlwiK2IudmFsdWUsYj1iW2NdLGImJihhKz1cIiA8IFwiKTtyZXR1cm4gYX19Zy5MUlVNYXA9bjtmdW5jdGlvbiBsKGEsYil7dGhpcy5rZXk9YSx0aGlzLnZhbHVlPWIsdGhpc1tjXT12b2lkIDAsdGhpc1tlXT12b2lkIDB9ZnVuY3Rpb24gaShhKXt0aGlzLmVudHJ5PWF9aS5wcm90b3R5cGVbU3ltYm9sLml0ZXJhdG9yXT1mdW5jdGlvbigpe3JldHVybiB0aGlzfSxpLnByb3RvdHlwZS5uZXh0PWZ1bmN0aW9uKCl7bGV0IGE9dGhpcy5lbnRyeTtyZXR1cm4gYT8odGhpcy5lbnRyeT1hW2NdLHtkb25lOiExLHZhbHVlOlthLmtleSxhLnZhbHVlXX0pOntkb25lOiEwLHZhbHVlOnZvaWQgMH19O2Z1bmN0aW9uIGooYSl7dGhpcy5lbnRyeT1hfWoucHJvdG90eXBlW1N5bWJvbC5pdGVyYXRvcl09ZnVuY3Rpb24oKXtyZXR1cm4gdGhpc30sai5wcm90b3R5cGUubmV4dD1mdW5jdGlvbigpe2xldCBhPXRoaXMuZW50cnk7cmV0dXJuIGE/KHRoaXMuZW50cnk9YVtjXSx7ZG9uZTohMSx2YWx1ZTphLmtleX0pOntkb25lOiEwLHZhbHVlOnZvaWQgMH19O2Z1bmN0aW9uIGsoYSl7dGhpcy5lbnRyeT1hfWsucHJvdG90eXBlW1N5bWJvbC5pdGVyYXRvcl09ZnVuY3Rpb24oKXtyZXR1cm4gdGhpc30say5wcm90b3R5cGUubmV4dD1mdW5jdGlvbigpe2xldCBhPXRoaXMuZW50cnk7cmV0dXJuIGE/KHRoaXMuZW50cnk9YVtjXSx7ZG9uZTohMSx2YWx1ZTphLnZhbHVlfSk6e2RvbmU6ITAsdmFsdWU6dm9pZCAwfX19KTtcbi8vIyBzb3VyY2VNYXBwaW5nVVJMPWxydS5qcy5tYXBcbiIsIlwidXNlIHN0cmljdFwiO1xuT2JqZWN0LmRlZmluZVByb3BlcnR5KGV4cG9ydHMsIFwiX19lc01vZHVsZVwiLCB7IHZhbHVlOiB0cnVlIH0pO1xuZXhwb3J0cy50aWxlTWFwT3B0aW9uc0NvbnRlbnRQcm9jZXNzb3IgPSBleHBvcnRzLlRpbGVNYXAgPSBleHBvcnRzLlRpbGVBbGlnbm1lbnQgPSB2b2lkIDA7XG5jb25zdCB2ZWNfMSA9IHJlcXVpcmUoXCJAYmFzZW1lbnR1bml2ZXJzZS92ZWNcIik7XG5jb25zdCBscnVfbWFwXzEgPSByZXF1aXJlKFwibHJ1X21hcFwiKTtcbnZhciBUaWxlQWxpZ25tZW50O1xuKGZ1bmN0aW9uIChUaWxlQWxpZ25tZW50KSB7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiVG9wTGVmdFwiXSA9IDBdID0gXCJUb3BMZWZ0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiVG9wXCJdID0gMV0gPSBcIlRvcFwiO1xuICAgIFRpbGVBbGlnbm1lbnRbVGlsZUFsaWdubWVudFtcIlRvcFJpZ2h0XCJdID0gMl0gPSBcIlRvcFJpZ2h0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiTGVmdFwiXSA9IDNdID0gXCJMZWZ0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiQ2VudGVyXCJdID0gNF0gPSBcIkNlbnRlclwiO1xuICAgIFRpbGVBbGlnbm1lbnRbVGlsZUFsaWdubWVudFtcIlJpZ2h0XCJdID0gNV0gPSBcIlJpZ2h0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiQm90dG9tTGVmdFwiXSA9IDZdID0gXCJCb3R0b21MZWZ0XCI7XG4gICAgVGlsZUFsaWdubWVudFtUaWxlQWxpZ25tZW50W1wiQm90dG9tXCJdID0gN10gPSBcIkJvdHRvbVwiO1xuICAgIFRpbGVBbGlnbm1lbnRbVGlsZUFsaWdubWVudFtcIkJvdHRvbVJpZ2h0XCJdID0gOF0gPSBcIkJvdHRvbVJpZ2h0XCI7XG59KShUaWxlQWxpZ25tZW50ID0gZXhwb3J0cy5UaWxlQWxpZ25tZW50IHx8IChleHBvcnRzLlRpbGVBbGlnbm1lbnQgPSB7fSkpO1xuZnVuY3Rpb24gY2xhbXAoYSwgbWluID0gMCwgbWF4ID0gMSkge1xuICAgIHJldHVybiBhIDwgbWluID8gbWluIDogKGEgPiBtYXggPyBtYXggOiBhKTtcbn1cbmZ1bmN0aW9uIHBvaW50SW5SZWN0YW5nbGUocG9pbnQsIHRvcExlZnQsIGJvdHRvbVJpZ2h0KSB7XG4gICAgcmV0dXJuIChwb2ludC54ID49IHRvcExlZnQueCAmJlxuICAgICAgICBwb2ludC55ID49IHRvcExlZnQueSAmJlxuICAgICAgICBwb2ludC54IDwgYm90dG9tUmlnaHQueCAmJlxuICAgICAgICBwb2ludC55IDwgYm90dG9tUmlnaHQueSk7XG59XG5jbGFzcyBUaWxlTWFwIHtcbiAgICBjb25zdHJ1Y3RvcihvcHRpb25zKSB7XG4gICAgICAgIGNvbnN0IGFjdHVhbE9wdGlvbnMgPSBPYmplY3QuYXNzaWduKHt9LCBUaWxlTWFwLkRFRkFVTFRfT1BUSU9OUywgb3B0aW9ucyAhPT0gbnVsbCAmJiBvcHRpb25zICE9PSB2b2lkIDAgPyBvcHRpb25zIDoge30pO1xuICAgICAgICBpZiAoIWFjdHVhbE9wdGlvbnMuZGVidWcgfHwgYWN0dWFsT3B0aW9ucy5kZWJ1ZyA9PT0gdHJ1ZSkge1xuICAgICAgICAgICAgYWN0dWFsT3B0aW9ucy5kZWJ1ZyA9IHtcbiAgICAgICAgICAgICAgICBzaG93T3JpZ2luOiAhIWFjdHVhbE9wdGlvbnMuZGVidWcsXG4gICAgICAgICAgICAgICAgc2hvd0NodW5rQm9yZGVyczogISFhY3R1YWxPcHRpb25zLmRlYnVnLFxuICAgICAgICAgICAgICAgIHNob3dDaHVua0xhYmVsczogISFhY3R1YWxPcHRpb25zLmRlYnVnLFxuICAgICAgICAgICAgICAgIHNob3dUaWxlQm9yZGVyczogISFhY3R1YWxPcHRpb25zLmRlYnVnLFxuICAgICAgICAgICAgfTtcbiAgICAgICAgfVxuICAgICAgICB0aGlzLm9wdGlvbnMgPSBhY3R1YWxPcHRpb25zO1xuICAgICAgICB0aGlzLmNodW5rQnVmZmVyID0gbmV3IGxydV9tYXBfMS5MUlVNYXAodGhpcy5vcHRpb25zLmNodW5rQnVmZmVyTWF4U2l6ZSk7XG4gICAgfVxuICAgIC8qKlxuICAgICAqIEdldCBhIG1pbmltYWwgc2V0IG9mIHJlY3RhbmdsZXMgd2hpY2ggY292ZXIgdGhlIHRpbGVzIGluIGEgZ2l2ZW4gbGF5ZXJcbiAgICAgKlxuICAgICAqIEBwYXJhbSBsYXllck5hbWUgVGhlIG5hbWUgb2YgdGhlIGxheWVyIHRvIGdldCByZWN0YW5nbGVzIGZvclxuICAgICAqIEBwYXJhbSBmaWVsZE5hbWUgV2Ugd2lsbCBjaGVjayB0aGUgdHJ1dGh5bmVzcyBvZiB0aGlzIGZpZWxkIGluIHRoZVxuICAgICAqIHRpbGUgZGVmaW5pdGlvblxuICAgICAqIEBwYXJhbSB0aWxlQm91bmRzIE9wdGlvbmFsIGJvdW5kcyB0byBjaGVja1xuICAgICAqL1xuICAgIGdldExheWVyUmVjdGFuZ2xlcyhsYXllck5hbWUsIGZpZWxkTmFtZSwgdGlsZUJvdW5kcykge1xuICAgICAgICAvLyBUT0RPXG4gICAgfVxuICAgIC8qKlxuICAgICAqIEdldCB0aGUgdGlsZSBhdCBhIGdpdmVuIHBvc2l0aW9uIGFuZCBpbiB0aGUgc3BlY2lmaWVkIGxheWVyXG4gICAgICpcbiAgICAgKiBJZiBubyBsYXllciBpcyBzcGVjaWZpZWQsIHJldHVybiBhIGRpY3Rpb25hcnkgb2YgbGF5ZXIgbmFtZXMgdG8gdGlsZVxuICAgICAqIGRlZmluaXRpb25zIChpLmUuIHJldHVybiBhbGwgbGF5ZXJzKVxuICAgICAqXG4gICAgICogSWYgbm8gdGlsZSBleGlzdHMgYXQgdGhpcyBwb3NpdGlvbiwgcmV0dXJuIG51bGxcbiAgICAgKi9cbiAgICBnZXRUaWxlQXRQb3NpdGlvbihwb3NpdGlvbiwgbGF5ZXJOYW1lKSB7XG4gICAgICAgIGlmIChsYXllck5hbWUpIHtcbiAgICAgICAgICAgIHJldHVybiB0aGlzLmdldFRpbGVBdFBvc2l0aW9uSW5MYXllcihwb3NpdGlvbiwgbGF5ZXJOYW1lKTtcbiAgICAgICAgfVxuICAgICAgICBjb25zdCByZXN1bHQgPSB7fTtcbiAgICAgICAgZm9yIChjb25zdCBsYXllciBvZiB0aGlzLm9wdGlvbnMubGF5ZXJzKSB7XG4gICAgICAgICAgICByZXN1bHRbbGF5ZXIubmFtZV0gPSB0aGlzLmdldFRpbGVBdFBvc2l0aW9uSW5MYXllcihwb3NpdGlvbiwgbGF5ZXIubmFtZSk7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHJlc3VsdDtcbiAgICB9XG4gICAgZ2V0VGlsZUF0UG9zaXRpb25JbkxheWVyKHBvc2l0aW9uLCBsYXllck5hbWUpIHtcbiAgICAgICAgdmFyIF9hLCBfYiwgX2M7XG4gICAgICAgIGNvbnN0IHRpbGVQb3NpdGlvbiA9IHZlY18xLnZlYy5tYXAodmVjXzEudmVjLm11bChwb3NpdGlvbiwgMSAvIHRoaXMub3B0aW9ucy50aWxlU2l6ZSksIE1hdGguZmxvb3IpO1xuICAgICAgICBjb25zdCBsYXllciA9IHRoaXMub3B0aW9ucy5sYXllcnMuZmluZCgobCkgPT4gbC5uYW1lID09PSBsYXllck5hbWUpO1xuICAgICAgICBpZiAoIWxheWVyKSB7XG4gICAgICAgICAgICByZXR1cm4gbnVsbDtcbiAgICAgICAgfVxuICAgICAgICBjb25zdCB0aWxlRGF0YSA9IChfYiA9IChfYSA9IGxheWVyLmRhdGEpID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYVt0aWxlUG9zaXRpb24ueV0pID09PSBudWxsIHx8IF9iID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYlt0aWxlUG9zaXRpb24ueF07XG4gICAgICAgIGlmICh0aWxlRGF0YSA9PT0gdW5kZWZpbmVkIHx8IHRpbGVEYXRhID09PSAtMSkge1xuICAgICAgICAgICAgcmV0dXJuIG51bGw7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKGxheWVyLnRpbGVzKSB7XG4gICAgICAgICAgICByZXR1cm4gKF9jID0gbGF5ZXIudGlsZXNbdGlsZURhdGFdKSAhPT0gbnVsbCAmJiBfYyAhPT0gdm9pZCAwID8gX2MgOiBudWxsO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBudWxsO1xuICAgIH1cbiAgICBoYXNoVmVjdG9yKHYpIHtcbiAgICAgICAgcmV0dXJuIHZlY18xLnZlYy5zdHIodik7XG4gICAgfVxuICAgIGRyYXcoY29udGV4dCwgc2NyZWVuLCBwb3NpdGlvbiwgc2NhbGUpIHtcbiAgICAgICAgdmFyIF9hLCBfYiwgX2MsIF9kO1xuICAgICAgICBjb25zdCBhYnNvbHV0ZUNodW5rU2l6ZSA9IHRoaXMub3B0aW9ucy50aWxlU2l6ZSAqIHRoaXMub3B0aW9ucy5jaHVua1NpemU7XG4gICAgICAgIGNvbnN0IGNodW5rQm9yZGVyID0gKDAsIHZlY18xLnZlYykodGhpcy5vcHRpb25zLmNodW5rQm9yZGVyKTtcbiAgICAgICAgLy8gTWF5YmUgY2xhbXAgc2NhbGVcbiAgICAgICAgbGV0IGFjdHVhbFNjYWxlID0gc2NhbGU7XG4gICAgICAgIGlmICh0aGlzLm9wdGlvbnMubWluU2NhbGUgJiYgYWN0dWFsU2NhbGUgPCB0aGlzLm9wdGlvbnMubWluU2NhbGUpIHtcbiAgICAgICAgICAgIGFjdHVhbFNjYWxlID0gdGhpcy5vcHRpb25zLm1pblNjYWxlO1xuICAgICAgICB9XG4gICAgICAgIGlmICh0aGlzLm9wdGlvbnMubWF4U2NhbGUgJiYgYWN0dWFsU2NhbGUgPiB0aGlzLm9wdGlvbnMubWF4U2NhbGUpIHtcbiAgICAgICAgICAgIGFjdHVhbFNjYWxlID0gdGhpcy5vcHRpb25zLm1heFNjYWxlO1xuICAgICAgICB9XG4gICAgICAgIC8vIE1heWJlIGNsYW1wIHBvc2l0aW9uIHRvIGJvdW5kc1xuICAgICAgICBsZXQgYWN0dWFsUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKShwb3NpdGlvbik7XG4gICAgICAgIGlmICh0aGlzLm9wdGlvbnMuYm91bmRzICYmIHRoaXMub3B0aW9ucy5jbGFtcFBvc2l0aW9uVG9Cb3VuZHMpIHtcbiAgICAgICAgICAgIGNvbnN0IHRpbGVTaXplU2NhbGVkID0gdGhpcy5vcHRpb25zLnRpbGVTaXplIC8gYWN0dWFsU2NhbGU7XG4gICAgICAgICAgICBjb25zdCBoYWxmU2NyZWVuU2NhbGVkID0gdmVjXzEudmVjLm1hcCh2ZWNfMS52ZWMubXVsKHNjcmVlbiwgMSAvIChhY3R1YWxTY2FsZSAqIDIpKSwgTWF0aC5jZWlsKTtcbiAgICAgICAgICAgIGNvbnN0IG1pblBvc2l0aW9uID0gKDAsIHZlY18xLnZlYykodGhpcy5vcHRpb25zLmJvdW5kcy50b3BMZWZ0LnggKiB0aWxlU2l6ZVNjYWxlZCArIGhhbGZTY3JlZW5TY2FsZWQueCwgdGhpcy5vcHRpb25zLmJvdW5kcy50b3BMZWZ0LnkgKiB0aWxlU2l6ZVNjYWxlZCArIGhhbGZTY3JlZW5TY2FsZWQueSk7XG4gICAgICAgICAgICBjb25zdCBtYXhQb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKHRoaXMub3B0aW9ucy5ib3VuZHMuYm90dG9tUmlnaHQueCAqIHRpbGVTaXplU2NhbGVkIC0gaGFsZlNjcmVlblNjYWxlZC54LCB0aGlzLm9wdGlvbnMuYm91bmRzLmJvdHRvbVJpZ2h0LnkgKiB0aWxlU2l6ZVNjYWxlZCAtIGhhbGZTY3JlZW5TY2FsZWQueSk7XG4gICAgICAgICAgICBhY3R1YWxQb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKGNsYW1wKGFjdHVhbFBvc2l0aW9uLngsIG1pblBvc2l0aW9uLngsIG1heFBvc2l0aW9uLngpLCBjbGFtcChhY3R1YWxQb3NpdGlvbi55LCBtaW5Qb3NpdGlvbi55LCBtYXhQb3NpdGlvbi55KSk7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3Qgc2NyZWVuU2l6ZUluQ2h1bmtzID0gdmVjXzEudmVjLm1hcCh2ZWNfMS52ZWMubXVsKHNjcmVlbiwgMSAvIChhYnNvbHV0ZUNodW5rU2l6ZSAqIGFjdHVhbFNjYWxlKSksIE1hdGguY2VpbCk7XG4gICAgICAgIGNvbnN0IHNjcmVlbkNlbnRlckNodW5rID0gdmVjXzEudmVjLm1hcCh2ZWNfMS52ZWMubXVsKGFjdHVhbFBvc2l0aW9uLCAxIC8gYWJzb2x1dGVDaHVua1NpemUpLCBNYXRoLmZsb29yKTtcbiAgICAgICAgY29uc3QgdG9wTGVmdENodW5rID0gdmVjXzEudmVjLnN1Yih2ZWNfMS52ZWMuc3ViKHNjcmVlbkNlbnRlckNodW5rLCB2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwoc2NyZWVuU2l6ZUluQ2h1bmtzLCAwLjUpLCBNYXRoLmNlaWwpKSwgY2h1bmtCb3JkZXIpO1xuICAgICAgICBjb25zdCBib3R0b21SaWdodENodW5rID0gdmVjXzEudmVjLmFkZCh2ZWNfMS52ZWMuYWRkKHNjcmVlbkNlbnRlckNodW5rLCB2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwoc2NyZWVuU2l6ZUluQ2h1bmtzLCAwLjUpLCBNYXRoLmNlaWwpKSwgY2h1bmtCb3JkZXIpO1xuICAgICAgICBjb250ZXh0LnNhdmUoKTtcbiAgICAgICAgY29udGV4dC5zY2FsZShhY3R1YWxTY2FsZSwgYWN0dWFsU2NhbGUpO1xuICAgICAgICBjb250ZXh0LnRyYW5zbGF0ZSgtYWN0dWFsUG9zaXRpb24ueCArIHNjcmVlbi54IC8gKGFjdHVhbFNjYWxlICogMiksIC1hY3R1YWxQb3NpdGlvbi55ICsgc2NyZWVuLnkgLyAoYWN0dWFsU2NhbGUgKiAyKSk7XG4gICAgICAgIChfYiA9IChfYSA9IHRoaXMub3B0aW9ucykucHJlUmVuZGVyKSA9PT0gbnVsbCB8fCBfYiA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2IuY2FsbChfYSwgY29udGV4dCwgdGhpcywgc2NyZWVuLCBhY3R1YWxQb3NpdGlvbiwgYWN0dWFsU2NhbGUpO1xuICAgICAgICAvLyBSZW5kZXIgY2h1bmtzXG4gICAgICAgIGZvciAobGV0IHkgPSB0b3BMZWZ0Q2h1bmsueTsgeSA8IGJvdHRvbVJpZ2h0Q2h1bmsueTsgeSsrKSB7XG4gICAgICAgICAgICBmb3IgKGxldCB4ID0gdG9wTGVmdENodW5rLng7IHggPCBib3R0b21SaWdodENodW5rLng7IHgrKykge1xuICAgICAgICAgICAgICAgIGNvbnN0IGNodW5rUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh4LCB5KTtcbiAgICAgICAgICAgICAgICBjb25zdCBjaHVua0Fic29sdXRlUG9zaXRpb24gPSB2ZWNfMS52ZWMubXVsKGNodW5rUG9zaXRpb24sIGFic29sdXRlQ2h1bmtTaXplKTtcbiAgICAgICAgICAgICAgICAvLyBDaGVjayBpZiB3ZSBoYXZlIHRoaXMgY2h1bmsgaW4gdGhlIGNhY2hlXG4gICAgICAgICAgICAgICAgY29uc3QgY2h1bmtIYXNoID0gdGhpcy5oYXNoVmVjdG9yKGNodW5rUG9zaXRpb24pO1xuICAgICAgICAgICAgICAgIGlmICghdGhpcy5jaHVua0J1ZmZlci5oYXMoY2h1bmtIYXNoKSkge1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmNodW5rQnVmZmVyLnNldChjaHVua0hhc2gsIHRoaXMuZ2VuZXJhdGVDaHVuayhjaHVua1Bvc2l0aW9uLCBhYnNvbHV0ZUNodW5rU2l6ZSkpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBjb25zdCBjaHVuayA9IHRoaXMuY2h1bmtCdWZmZXIuZ2V0KGNodW5rSGFzaCk7XG4gICAgICAgICAgICAgICAgaWYgKGNodW5rKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnRleHQuZHJhd0ltYWdlKGNodW5rLmltYWdlLCBjaHVua0Fic29sdXRlUG9zaXRpb24ueCwgY2h1bmtBYnNvbHV0ZVBvc2l0aW9uLnkpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICAoX2QgPSAoX2MgPSB0aGlzLm9wdGlvbnMpLnBvc3RSZW5kZXIpID09PSBudWxsIHx8IF9kID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfZC5jYWxsKF9jLCBjb250ZXh0LCB0aGlzLCBzY3JlZW4sIGFjdHVhbFBvc2l0aW9uLCBhY3R1YWxTY2FsZSk7XG4gICAgICAgIC8vIFJlbmRlciBkZWJ1ZyBoZWxwZXJzXG4gICAgICAgIGlmICh0aGlzLm9wdGlvbnMuZGVidWcuc2hvd1RpbGVCb3JkZXJzKSB7XG4gICAgICAgICAgICBjb25zdCB0b3BMZWZ0VGlsZSA9IHZlY18xLnZlYy5tdWwodmVjXzEudmVjLnN1YihzY3JlZW5DZW50ZXJDaHVuaywgdmVjXzEudmVjLmFkZCh2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwoc2NyZWVuU2l6ZUluQ2h1bmtzLCAwLjUpLCBNYXRoLmNlaWwpLCAoMCwgdmVjXzEudmVjKSgxKSkpLCB0aGlzLm9wdGlvbnMuY2h1bmtTaXplKTtcbiAgICAgICAgICAgIGNvbnN0IGJvdHRvbVJpZ2h0VGlsZSA9IHZlY18xLnZlYy5tdWwodmVjXzEudmVjLmFkZChzY3JlZW5DZW50ZXJDaHVuaywgdmVjXzEudmVjLmFkZCh2ZWNfMS52ZWMubWFwKHZlY18xLnZlYy5tdWwoc2NyZWVuU2l6ZUluQ2h1bmtzLCAwLjUpLCBNYXRoLmNlaWwpLCAoMCwgdmVjXzEudmVjKSgxKSkpLCB0aGlzLm9wdGlvbnMuY2h1bmtTaXplKTtcbiAgICAgICAgICAgIGZvciAobGV0IHkgPSB0b3BMZWZ0VGlsZS55OyB5IDwgYm90dG9tUmlnaHRUaWxlLnk7IHkrKykge1xuICAgICAgICAgICAgICAgIHRoaXMuZHJhd0xpbmUoY29udGV4dCwgKDAsIHZlY18xLnZlYykoYWN0dWFsUG9zaXRpb24ueCAtIHNjcmVlbi54IC8gKGFjdHVhbFNjYWxlICogMiksIHkgKiB0aGlzLm9wdGlvbnMudGlsZVNpemUpLCAoMCwgdmVjXzEudmVjKShhY3R1YWxQb3NpdGlvbi54ICsgc2NyZWVuLnggLyAoYWN0dWFsU2NhbGUgKiAyKSwgeSAqIHRoaXMub3B0aW9ucy50aWxlU2l6ZSksIFRpbGVNYXAuREVCVUdfVElMRV9CT1JERVJfQ09MT1VSLCBUaWxlTWFwLkRFQlVHX1RJTEVfQk9SREVSX0xJTkVfV0lEVEgpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZm9yIChsZXQgeCA9IHRvcExlZnRUaWxlLng7IHggPCBib3R0b21SaWdodFRpbGUueDsgeCsrKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5kcmF3TGluZShjb250ZXh0LCAoMCwgdmVjXzEudmVjKSh4ICogdGhpcy5vcHRpb25zLnRpbGVTaXplLCBhY3R1YWxQb3NpdGlvbi55IC0gc2NyZWVuLnkgLyAoYWN0dWFsU2NhbGUgKiAyKSksICgwLCB2ZWNfMS52ZWMpKHggKiB0aGlzLm9wdGlvbnMudGlsZVNpemUsIGFjdHVhbFBvc2l0aW9uLnkgKyBzY3JlZW4ueSAvIChhY3R1YWxTY2FsZSAqIDIpKSwgVGlsZU1hcC5ERUJVR19USUxFX0JPUkRFUl9DT0xPVVIsIFRpbGVNYXAuREVCVUdfVElMRV9CT1JERVJfTElORV9XSURUSCk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5kZWJ1Zy5zaG93Q2h1bmtCb3JkZXJzKSB7XG4gICAgICAgICAgICBmb3IgKGxldCB5ID0gdG9wTGVmdENodW5rLnk7IHkgPCBib3R0b21SaWdodENodW5rLnk7IHkrKykge1xuICAgICAgICAgICAgICAgIHRoaXMuZHJhd0xpbmUoY29udGV4dCwgKDAsIHZlY18xLnZlYykoYWN0dWFsUG9zaXRpb24ueCAtIHNjcmVlbi54IC8gKGFjdHVhbFNjYWxlICogMiksIHkgKiBhYnNvbHV0ZUNodW5rU2l6ZSksICgwLCB2ZWNfMS52ZWMpKGFjdHVhbFBvc2l0aW9uLnggKyBzY3JlZW4ueCAvIChhY3R1YWxTY2FsZSAqIDIpLCB5ICogYWJzb2x1dGVDaHVua1NpemUpLCBUaWxlTWFwLkRFQlVHX0NIVU5LX0JPUkRFUl9DT0xPVVIsIFRpbGVNYXAuREVCVUdfQ0hVTktfQk9SREVSX0xJTkVfV0lEVEgpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZm9yIChsZXQgeCA9IHRvcExlZnRDaHVuay54OyB4IDwgYm90dG9tUmlnaHRDaHVuay54OyB4KyspIHtcbiAgICAgICAgICAgICAgICB0aGlzLmRyYXdMaW5lKGNvbnRleHQsICgwLCB2ZWNfMS52ZWMpKHggKiBhYnNvbHV0ZUNodW5rU2l6ZSwgYWN0dWFsUG9zaXRpb24ueSAtIHNjcmVlbi55IC8gKGFjdHVhbFNjYWxlICogMikpLCAoMCwgdmVjXzEudmVjKSh4ICogYWJzb2x1dGVDaHVua1NpemUsIGFjdHVhbFBvc2l0aW9uLnkgKyBzY3JlZW4ueSAvIChhY3R1YWxTY2FsZSAqIDIpKSwgVGlsZU1hcC5ERUJVR19DSFVOS19CT1JERVJfQ09MT1VSLCBUaWxlTWFwLkRFQlVHX0NIVU5LX0JPUkRFUl9MSU5FX1dJRFRIKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBpZiAodGhpcy5vcHRpb25zLmRlYnVnLnNob3dDaHVua0xhYmVscykge1xuICAgICAgICAgICAgY29udGV4dC5zYXZlKCk7XG4gICAgICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFRpbGVNYXAuREVCVUdfQ0hVTktfTEFCRUxfQ09MT1VSO1xuICAgICAgICAgICAgY29udGV4dC5mb250ID0gVGlsZU1hcC5ERUJVR19DSFVOS19MQUJFTF9GT05UO1xuICAgICAgICAgICAgY29udGV4dC50ZXh0QmFzZWxpbmUgPSAnbWlkZGxlJztcbiAgICAgICAgICAgIGNvbnRleHQudGV4dEFsaWduID0gJ2NlbnRlcic7XG4gICAgICAgICAgICBmb3IgKGxldCB5ID0gdG9wTGVmdENodW5rLnk7IHkgPCBib3R0b21SaWdodENodW5rLnk7IHkrKykge1xuICAgICAgICAgICAgICAgIGZvciAobGV0IHggPSB0b3BMZWZ0Q2h1bmsueDsgeCA8IGJvdHRvbVJpZ2h0Q2h1bmsueDsgeCsrKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnRleHQuZmlsbFRleHQoYCR7eH0sICR7eX1gLCB4ICogYWJzb2x1dGVDaHVua1NpemUgKyBhYnNvbHV0ZUNodW5rU2l6ZSAvIDIsIHkgKiBhYnNvbHV0ZUNodW5rU2l6ZSArIGFic29sdXRlQ2h1bmtTaXplIC8gMik7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgY29udGV4dC5yZXN0b3JlKCk7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5kZWJ1Zy5zaG93T3JpZ2luICYmXG4gICAgICAgICAgICBwb2ludEluUmVjdGFuZ2xlKCgwLCB2ZWNfMS52ZWMpKDAsIDApLCB0b3BMZWZ0Q2h1bmssIGJvdHRvbVJpZ2h0Q2h1bmspKSB7XG4gICAgICAgICAgICB0aGlzLmRyYXdDcm9zcyhjb250ZXh0LCAoMCwgdmVjXzEudmVjKSgwLCAwKSwgVGlsZU1hcC5ERUJVR19PUklHSU5fQ09MT1VSLCBUaWxlTWFwLkRFQlVHX09SSUdJTl9MSU5FX1dJRFRILCBUaWxlTWFwLkRFQlVHX09SSUdJTl9TSVpFKTtcbiAgICAgICAgfVxuICAgICAgICBjb250ZXh0LnJlc3RvcmUoKTtcbiAgICB9XG4gICAgZ2VuZXJhdGVDaHVuayhjaHVua1Bvc2l0aW9uLCBhYnNvbHV0ZUNodW5rU2l6ZSkge1xuICAgICAgICB2YXIgX2EsIF9iLCBfYywgX2QsIF9lLCBfZiwgX2csIF9oLCBfaiwgX2ssIF9sLCBfbTtcbiAgICAgICAgY29uc3QgY2h1bmtDYW52YXMgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KCdjYW52YXMnKTtcbiAgICAgICAgY29uc3QgY2h1bmtDb250ZXh0ID0gY2h1bmtDYW52YXMuZ2V0Q29udGV4dCgnMmQnKTtcbiAgICAgICAgY2h1bmtDYW52YXMud2lkdGggPSBhYnNvbHV0ZUNodW5rU2l6ZTtcbiAgICAgICAgY2h1bmtDYW52YXMuaGVpZ2h0ID0gYWJzb2x1dGVDaHVua1NpemU7XG4gICAgICAgIGxldCBjaHVuayA9IHtcbiAgICAgICAgICAgIGNodW5rUG9zaXRpb24sXG4gICAgICAgICAgICBpbWFnZTogY2h1bmtDYW52YXMsXG4gICAgICAgIH07XG4gICAgICAgIGNvbnN0IHRvcExlZnRUaWxlID0gdmVjXzEudmVjLm11bChjaHVua1Bvc2l0aW9uLCB0aGlzLm9wdGlvbnMuY2h1bmtTaXplKTtcbiAgICAgICAgY29uc3QgYm90dG9tUmlnaHRUaWxlID0gdmVjXzEudmVjLmFkZCh0b3BMZWZ0VGlsZSwgKDAsIHZlY18xLnZlYykodGhpcy5vcHRpb25zLmNodW5rU2l6ZSAtIDEpKTtcbiAgICAgICAgY29uc3QgYm91bmRzVG9wTGVmdCA9IChfYiA9IChfYSA9IHRoaXMub3B0aW9ucy5ib3VuZHMpID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS50b3BMZWZ0KSAhPT0gbnVsbCAmJiBfYiAhPT0gdm9pZCAwID8gX2IgOiAoMCwgdmVjXzEudmVjKSgwKTtcbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5wcmVHZW5lcmF0ZUNodW5rKSB7XG4gICAgICAgICAgICBjb25zdCByZXN1bHQgPSB0aGlzLm9wdGlvbnMucHJlR2VuZXJhdGVDaHVuayhjaHVua0NvbnRleHQsIHRoaXMsIHtcbiAgICAgICAgICAgICAgICB0b3BMZWZ0OiB0b3BMZWZ0VGlsZSxcbiAgICAgICAgICAgICAgICBib3R0b21SaWdodDogYm90dG9tUmlnaHRUaWxlLFxuICAgICAgICAgICAgfSwgY2h1bmtQb3NpdGlvbik7XG4gICAgICAgICAgICBpZiAoQXJyYXkuaXNBcnJheShyZXN1bHQpKSB7XG4gICAgICAgICAgICAgICAgaWYgKCFyZXN1bHRbMV0pIHtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNodW5rO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICAvLyBEZWZhdWx0IGdlbmVyYXRpb24sIHJlbmRlciB0aWxlcyBmcm9tIHRpbGVtYXAgZGF0YVxuICAgICAgICBmb3IgKGNvbnN0IGxheWVyIG9mIHRoaXMub3B0aW9ucy5sYXllcnMpIHtcbiAgICAgICAgICAgIGNodW5rQ29udGV4dC5zYXZlKCk7XG4gICAgICAgICAgICBjaHVua0NvbnRleHQuZ2xvYmFsQWxwaGEgPSAoX2MgPSBsYXllci5vcGFjaXR5KSAhPT0gbnVsbCAmJiBfYyAhPT0gdm9pZCAwID8gX2MgOiAxO1xuICAgICAgICAgICAgY29uc3QgYWxpZ25tZW50ID0gKF9kID0gbGF5ZXIuYWxpZ25tZW50KSAhPT0gbnVsbCAmJiBfZCAhPT0gdm9pZCAwID8gX2QgOiBUaWxlQWxpZ25tZW50LkNlbnRlcjtcbiAgICAgICAgICAgIGZvciAobGV0IHkgPSB0b3BMZWZ0VGlsZS55OyB5IDw9IGJvdHRvbVJpZ2h0VGlsZS55OyB5KyspIHtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCB4ID0gdG9wTGVmdFRpbGUueDsgeCA8PSBib3R0b21SaWdodFRpbGUueDsgeCsrKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IHRpbGVQb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKHgsIHkpO1xuICAgICAgICAgICAgICAgICAgICAoX2UgPSBsYXllci5wcmVSZW5kZXJUaWxlKSA9PT0gbnVsbCB8fCBfZSA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2UuY2FsbChsYXllciwgY2h1bmtDb250ZXh0LCB0aGlzLCBsYXllciwgY2h1bmtQb3NpdGlvbiwgdGlsZVBvc2l0aW9uKTtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgdGlsZURhdGFQb3NpdGlvbiA9IHZlY18xLnZlYy5zdWIodGlsZVBvc2l0aW9uLCBib3VuZHNUb3BMZWZ0KTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKHRpbGVEYXRhUG9zaXRpb24ueCA8IDAgfHwgdGlsZURhdGFQb3NpdGlvbi55IDwgMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgY29uc3QgdGlsZURhdGEgPSAoX2cgPSAoX2YgPSBsYXllci5kYXRhKSA9PT0gbnVsbCB8fCBfZiA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2ZbdGlsZURhdGFQb3NpdGlvbi55XSkgPT09IG51bGwgfHwgX2cgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9nW3RpbGVEYXRhUG9zaXRpb24ueF07XG4gICAgICAgICAgICAgICAgICAgIGlmICh0aWxlRGF0YSA9PT0gdW5kZWZpbmVkIHx8IHRpbGVEYXRhID09PSAtMSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgY29uc3QgdGlsZUltYWdlID0gKF9qID0gKF9oID0gbGF5ZXIudGlsZXMpID09PSBudWxsIHx8IF9oID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfaFt0aWxlRGF0YV0pID09PSBudWxsIHx8IF9qID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfai5pbWFnZTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKCF0aWxlSW1hZ2UpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IHRpbGVBYnNvbHV0ZVBvc2l0aW9uID0gdmVjXzEudmVjLnN1Yih2ZWNfMS52ZWMubXVsKHRpbGVQb3NpdGlvbiwgdGhpcy5vcHRpb25zLnRpbGVTaXplKSwgdmVjXzEudmVjLm11bChjaHVua1Bvc2l0aW9uLCBhYnNvbHV0ZUNodW5rU2l6ZSkpO1xuICAgICAgICAgICAgICAgICAgICAvLyBUaWxlIGNsaXBwaW5nXG4gICAgICAgICAgICAgICAgICAgIGlmIChsYXllci5jbGlwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBjaHVua0NvbnRleHQuc2F2ZSgpO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2h1bmtDb250ZXh0LmJlZ2luUGF0aCgpO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2h1bmtDb250ZXh0LnJlY3QodGlsZUFic29sdXRlUG9zaXRpb24ueCwgdGlsZUFic29sdXRlUG9zaXRpb24ueSwgdGhpcy5vcHRpb25zLnRpbGVTaXplLCB0aGlzLm9wdGlvbnMudGlsZVNpemUpO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2h1bmtDb250ZXh0LmNsaXAoKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAvLyBUaWxlIGFsaWdubWVudFxuICAgICAgICAgICAgICAgICAgICBsZXQgdGlsZUltYWdlQWJzb2x1dGVQb3NpdGlvbjtcbiAgICAgICAgICAgICAgICAgICAgc3dpdGNoIChhbGlnbm1lbnQpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNhc2UgVGlsZUFsaWdubWVudC5Ub3BMZWZ0OlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh0aWxlQWJzb2x1dGVQb3NpdGlvbik7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFRpbGVBbGlnbm1lbnQuVG9wOlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSgodGlsZUFic29sdXRlUG9zaXRpb24ueCArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAvIDIpIC0gdGlsZUltYWdlLndpZHRoIC8gMiwgdGlsZUFic29sdXRlUG9zaXRpb24ueSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFRpbGVBbGlnbm1lbnQuVG9wUmlnaHQ6XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdGlsZUltYWdlQWJzb2x1dGVQb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKHRpbGVBYnNvbHV0ZVBvc2l0aW9uLnggKyB0aGlzLm9wdGlvbnMudGlsZVNpemUgLSB0aWxlSW1hZ2Uud2lkdGgsIHRpbGVBYnNvbHV0ZVBvc2l0aW9uLnkpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2FzZSBUaWxlQWxpZ25tZW50LkxlZnQ6XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdGlsZUltYWdlQWJzb2x1dGVQb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKHRpbGVBYnNvbHV0ZVBvc2l0aW9uLngsICh0aWxlQWJzb2x1dGVQb3NpdGlvbi55ICsgdGhpcy5vcHRpb25zLnRpbGVTaXplIC8gMikgLSB0aWxlSW1hZ2UuaGVpZ2h0IC8gMik7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFRpbGVBbGlnbm1lbnQuQ2VudGVyOlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSgodGlsZUFic29sdXRlUG9zaXRpb24ueCArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAvIDIpIC0gdGlsZUltYWdlLndpZHRoIC8gMiwgKHRpbGVBYnNvbHV0ZVBvc2l0aW9uLnkgKyB0aGlzLm9wdGlvbnMudGlsZVNpemUgLyAyKSAtIHRpbGVJbWFnZS5oZWlnaHQgLyAyKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICAgICAgICAgIGNhc2UgVGlsZUFsaWdubWVudC5SaWdodDpcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB0aWxlSW1hZ2VBYnNvbHV0ZVBvc2l0aW9uID0gKDAsIHZlY18xLnZlYykodGlsZUFic29sdXRlUG9zaXRpb24ueCArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAtIHRpbGVJbWFnZS53aWR0aCwgKHRpbGVBYnNvbHV0ZVBvc2l0aW9uLnkgKyB0aGlzLm9wdGlvbnMudGlsZVNpemUgLyAyKSAtIHRpbGVJbWFnZS5oZWlnaHQgLyAyKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICAgICAgICAgIGNhc2UgVGlsZUFsaWdubWVudC5Cb3R0b21MZWZ0OlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh0aWxlQWJzb2x1dGVQb3NpdGlvbi54LCB0aWxlQWJzb2x1dGVQb3NpdGlvbi55ICsgdGhpcy5vcHRpb25zLnRpbGVTaXplIC0gdGlsZUltYWdlLmhlaWdodCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFRpbGVBbGlnbm1lbnQuQm90dG9tOlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSgodGlsZUFic29sdXRlUG9zaXRpb24ueCArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAvIDIpIC0gdGlsZUltYWdlLndpZHRoIC8gMiwgdGlsZUFic29sdXRlUG9zaXRpb24ueSArIHRoaXMub3B0aW9ucy50aWxlU2l6ZSAtIHRpbGVJbWFnZS5oZWlnaHQpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2FzZSBUaWxlQWxpZ25tZW50LkJvdHRvbVJpZ2h0OlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24gPSAoMCwgdmVjXzEudmVjKSh0aWxlQWJzb2x1dGVQb3NpdGlvbi54ICsgdGhpcy5vcHRpb25zLnRpbGVTaXplIC0gdGlsZUltYWdlLndpZHRoLCB0aWxlQWJzb2x1dGVQb3NpdGlvbi55ICsgdGhpcy5vcHRpb25zLnRpbGVTaXplIC0gdGlsZUltYWdlLmhlaWdodCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgY2h1bmtDb250ZXh0LmRyYXdJbWFnZSh0aWxlSW1hZ2UsIHRpbGVJbWFnZUFic29sdXRlUG9zaXRpb24ueCwgdGlsZUltYWdlQWJzb2x1dGVQb3NpdGlvbi55KTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGxheWVyLmNsaXApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNodW5rQ29udGV4dC5yZXN0b3JlKCk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgKF9rID0gbGF5ZXIucG9zdFJlbmRlclRpbGUpID09PSBudWxsIHx8IF9rID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfay5jYWxsKGxheWVyLCBjaHVua0NhbnZhcywgY2h1bmtDb250ZXh0LCB0aGlzLCBsYXllciwgY2h1bmtQb3NpdGlvbiwgdGlsZVBvc2l0aW9uKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBjaHVua0NvbnRleHQucmVzdG9yZSgpO1xuICAgICAgICB9XG4gICAgICAgIChfbSA9IChfbCA9IHRoaXMub3B0aW9ucykucG9zdEdlbmVyYXRlQ2h1bmspID09PSBudWxsIHx8IF9tID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfbS5jYWxsKF9sLCBjaHVua0NhbnZhcywgY2h1bmtDb250ZXh0LCB0aGlzLCB7XG4gICAgICAgICAgICB0b3BMZWZ0OiB0b3BMZWZ0VGlsZSxcbiAgICAgICAgICAgIGJvdHRvbVJpZ2h0OiBib3R0b21SaWdodFRpbGUsXG4gICAgICAgIH0sIGNodW5rUG9zaXRpb24pO1xuICAgICAgICByZXR1cm4gY2h1bms7XG4gICAgfVxuICAgIGRyYXdMaW5lKGNvbnRleHQsIHN0YXJ0LCBlbmQsIGNvbG91ciwgbGluZVdpZHRoKSB7XG4gICAgICAgIGNvbnRleHQuc2F2ZSgpO1xuICAgICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IGxpbmVXaWR0aDtcbiAgICAgICAgY29udGV4dC5zdHJva2VTdHlsZSA9IGNvbG91cjtcbiAgICAgICAgY29udGV4dC5iZWdpblBhdGgoKTtcbiAgICAgICAgY29udGV4dC5tb3ZlVG8oc3RhcnQueCwgc3RhcnQueSk7XG4gICAgICAgIGNvbnRleHQubGluZVRvKGVuZC54LCBlbmQueSk7XG4gICAgICAgIGNvbnRleHQuc3Ryb2tlKCk7XG4gICAgICAgIGNvbnRleHQucmVzdG9yZSgpO1xuICAgIH1cbiAgICBkcmF3Q3Jvc3MoY29udGV4dCwgcG9zaXRpb24sIGNvbG91ciwgbGluZVdpZHRoLCBzaXplKSB7XG4gICAgICAgIGNvbnRleHQuc2F2ZSgpO1xuICAgICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IGxpbmVXaWR0aDtcbiAgICAgICAgY29uc3QgaGFsZlNpemUgPSBNYXRoLmNlaWwoc2l6ZSAvIDIpO1xuICAgICAgICBjb250ZXh0LnN0cm9rZVN0eWxlID0gY29sb3VyO1xuICAgICAgICBjb250ZXh0LmJlZ2luUGF0aCgpO1xuICAgICAgICBjb250ZXh0Lm1vdmVUbyhwb3NpdGlvbi54IC0gaGFsZlNpemUsIHBvc2l0aW9uLnkgLSBoYWxmU2l6ZSk7XG4gICAgICAgIGNvbnRleHQubGluZVRvKHBvc2l0aW9uLnggKyBoYWxmU2l6ZSwgcG9zaXRpb24ueSArIGhhbGZTaXplKTtcbiAgICAgICAgY29udGV4dC5tb3ZlVG8ocG9zaXRpb24ueCAtIGhhbGZTaXplLCBwb3NpdGlvbi55ICsgaGFsZlNpemUpO1xuICAgICAgICBjb250ZXh0LmxpbmVUbyhwb3NpdGlvbi54ICsgaGFsZlNpemUsIHBvc2l0aW9uLnkgLSBoYWxmU2l6ZSk7XG4gICAgICAgIGNvbnRleHQuc3Ryb2tlKCk7XG4gICAgICAgIGNvbnRleHQucmVzdG9yZSgpO1xuICAgIH1cbn1cbmV4cG9ydHMuVGlsZU1hcCA9IFRpbGVNYXA7XG4vLyBUT0RPXG4vLyBodHRwczovL3d3dy5ucG1qcy5jb20vcGFja2FnZS9yZWN0YW5nbGUtZGVjb21wb3NpdGlvblxuLy8gaHR0cHM6Ly93d3cubnBtanMuY29tL3BhY2thZ2UvZmFzdC1ybGVcblRpbGVNYXAuREVGQVVMVF9PUFRJT05TID0ge1xuICAgIGNsYW1wUG9zaXRpb25Ub0JvdW5kczogdHJ1ZSxcbiAgICB0aWxlU2l6ZTogMTYsXG4gICAgbGF5ZXJzOiBbXG4gICAgICAgIHtcbiAgICAgICAgICAgIG5hbWU6ICdkZWZhdWx0JyxcbiAgICAgICAgfSxcbiAgICBdLFxuICAgIGNodW5rU2l6ZTogOCxcbiAgICBjaHVua0JvcmRlcjogMSxcbiAgICBjaHVua0J1ZmZlck1heFNpemU6IDY0LFxufTtcblRpbGVNYXAuREVCVUdfT1JJR0lOX0NPTE9VUiA9ICdjeWFuJztcblRpbGVNYXAuREVCVUdfT1JJR0lOX0xJTkVfV0lEVEggPSAyO1xuVGlsZU1hcC5ERUJVR19PUklHSU5fU0laRSA9IDEwO1xuVGlsZU1hcC5ERUJVR19DSFVOS19CT1JERVJfQ09MT1VSID0gJ3llbGxvdyc7XG5UaWxlTWFwLkRFQlVHX0NIVU5LX0JPUkRFUl9MSU5FX1dJRFRIID0gMjtcblRpbGVNYXAuREVCVUdfQ0hVTktfTEFCRUxfQ09MT1VSID0gJ3doaXRlJztcblRpbGVNYXAuREVCVUdfQ0hVTktfTEFCRUxfRk9OVCA9ICcxMnB4IG1vbm9zcGFjZSc7XG5UaWxlTWFwLkRFQlVHX1RJTEVfQk9SREVSX0NPTE9VUiA9ICdvcmFuZ2UnO1xuVGlsZU1hcC5ERUJVR19USUxFX0JPUkRFUl9MSU5FX1dJRFRIID0gMTtcbi8qKlxuICogQ29udGVudCBNYW5hZ2VyIFByb2Nlc3NvciB3cmFwcGVyIHdoaWNoIGNvbnZlcnRzIFRpbGVNYXBPcHRpb25zRGF0YSBpbnRvXG4gKiBUaWxlTWFwT3B0aW9uc1xuICpcbiAqIEBzZWUgaHR0cHM6Ly93d3cubnBtanMuY29tL3BhY2thZ2UvQGJhc2VtZW50dW5pdmVyc2UvY29udGVudC1tYW5hZ2VyXG4gKi9cbmFzeW5jIGZ1bmN0aW9uIHRpbGVNYXBPcHRpb25zQ29udGVudFByb2Nlc3Nvcihjb250ZW50LCBkYXRhKSB7XG4gICAgLy9cbn1cbmV4cG9ydHMudGlsZU1hcE9wdGlvbnNDb250ZW50UHJvY2Vzc29yID0gdGlsZU1hcE9wdGlvbnNDb250ZW50UHJvY2Vzc29yO1xuLy8jIHNvdXJjZU1hcHBpbmdVUkw9ZGF0YTphcHBsaWNhdGlvbi9qc29uO2Jhc2U2NCxleUoyWlhKemFXOXVJam96TENKbWFXeGxJam9pYVc1a1pYZ3Vhbk1pTENKemIzVnlZMlZTYjI5MElqb2lJaXdpYzI5MWNtTmxjeUk2V3lJdUxpOXBibVJsZUM1MGN5SmRMQ0p1WVcxbGN5STZXMTBzSW0xaGNIQnBibWR6SWpvaU96czdRVUZCUVN3clEwRkJORU03UVVGRE5VTXNjVU5CUVdsRE8wRkJObEZxUXl4SlFVRlpMR0ZCVlZnN1FVRldSQ3hYUVVGWkxHRkJRV0U3U1VGRGRrSXNkVVJCUVZjc1EwRkJRVHRKUVVOWUxDdERRVUZITEVOQlFVRTdTVUZEU0N4NVJFRkJVU3hEUVVGQk8wbEJRMUlzYVVSQlFVa3NRMEZCUVR0SlFVTktMSEZFUVVGTkxFTkJRVUU3U1VGRFRpeHRSRUZCU3l4RFFVRkJPMGxCUTB3c05rUkJRVlVzUTBGQlFUdEpRVU5XTEhGRVFVRk5MRU5CUVVFN1NVRkRUaXdyUkVGQlZ5eERRVUZCTzBGQlEySXNRMEZCUXl4RlFWWlhMR0ZCUVdFc1IwRkJZaXh4UWtGQllTeExRVUZpTEhGQ1FVRmhMRkZCVlhoQ08wRkJZMFFzVTBGQlV5eExRVUZMTEVOQlFVTXNRMEZCVXl4RlFVRkZMRTFCUVdNc1EwRkJReXhGUVVGRkxFMUJRV01zUTBGQlF6dEpRVU40UkN4UFFVRlBMRU5CUVVNc1IwRkJSeXhIUVVGSExFTkJRVU1zUTBGQlF5eERRVUZETEVkQlFVY3NRMEZCUXl4RFFVRkRMRU5CUVVNc1EwRkJReXhEUVVGRExFZEJRVWNzUjBGQlJ5eERRVUZETEVOQlFVTXNRMEZCUXl4SFFVRkhMRU5CUVVNc1EwRkJReXhEUVVGRExFTkJRVU1zUTBGQlF5eERRVUZETzBGQlF6ZERMRU5CUVVNN1FVRkZSQ3hUUVVGVExHZENRVUZuUWl4RFFVTjJRaXhMUVVGVkxFVkJRMVlzVDBGQldTeEZRVU5hTEZkQlFXZENPMGxCUldoQ0xFOUJRVThzUTBGRFRDeExRVUZMTEVOQlFVTXNRMEZCUXl4SlFVRkpMRTlCUVU4c1EwRkJReXhEUVVGRE8xRkJRM0JDTEV0QlFVc3NRMEZCUXl4RFFVRkRMRWxCUVVrc1QwRkJUeXhEUVVGRExFTkJRVU03VVVGRGNFSXNTMEZCU3l4RFFVRkRMRU5CUVVNc1IwRkJSeXhYUVVGWExFTkJRVU1zUTBGQlF6dFJRVU4yUWl4TFFVRkxMRU5CUVVNc1EwRkJReXhIUVVGSExGZEJRVmNzUTBGQlF5eERRVUZETEVOQlEzaENMRU5CUVVNN1FVRkRTaXhEUVVGRE8wRkJSVVFzVFVGQllTeFBRVUZQTzBsQmMwTnNRaXhaUVVGdFFpeFBRVUZ2UXp0UlFVTnlSQ3hOUVVGTkxHRkJRV0VzUjBGQlJ5eE5RVUZOTEVOQlFVTXNUVUZCVFN4RFFVTnFReXhGUVVGRkxFVkJRMFlzVDBGQlR5eERRVUZETEdWQlFXVXNSVUZEZGtJc1QwRkJUeXhoUVVGUUxFOUJRVThzWTBGQlVDeFBRVUZQTEVkQlFVa3NSVUZCUlN4RFFVTmtMRU5CUVVNN1VVRkZSaXhKUVVGSkxFTkJRVU1zWVVGQllTeERRVUZETEV0QlFVc3NTVUZCU1N4aFFVRmhMRU5CUVVNc1MwRkJTeXhMUVVGTExFbEJRVWtzUlVGQlJUdFpRVU40UkN4aFFVRmhMRU5CUVVNc1MwRkJTeXhIUVVGSE8yZENRVU53UWl4VlFVRlZMRVZCUVVVc1EwRkJReXhEUVVGRExHRkJRV0VzUTBGQlF5eExRVUZMTzJkQ1FVTnFReXhuUWtGQlowSXNSVUZCUlN4RFFVRkRMRU5CUVVNc1lVRkJZU3hEUVVGRExFdEJRVXM3WjBKQlEzWkRMR1ZCUVdVc1JVRkJSU3hEUVVGRExFTkJRVU1zWVVGQllTeERRVUZETEV0QlFVczdaMEpCUTNSRExHVkJRV1VzUlVGQlJTeERRVUZETEVOQlFVTXNZVUZCWVN4RFFVRkRMRXRCUVVzN1lVRkRka01zUTBGQlF6dFRRVU5JTzFGQlJVUXNTVUZCU1N4RFFVRkRMRTlCUVU4c1IwRkJSeXhoUVVGdlF5eERRVUZETzFGQlJYQkVMRWxCUVVrc1EwRkJReXhYUVVGWExFZEJRVWNzU1VGQlNTeG5Ra0ZCVFN4RFFVRkRMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zYTBKQlFXdENMRU5CUVVNc1EwRkJRenRKUVVOcVJTeERRVUZETzBsQlJVUTdPenM3T3pzN1QwRlBSenRKUVVOSkxHdENRVUZyUWl4RFFVTjJRaXhUUVVGcFFpeEZRVU5xUWl4VFFVRnJRaXhGUVVOc1FpeFZRVUZ0UWp0UlFVVnVRaXhQUVVGUE8wbEJRMVFzUTBGQlF6dEpRVVZFT3pzN096czdPMDlCVDBjN1NVRkRTU3hwUWtGQmFVSXNRMEZEZEVJc1VVRkJZU3hGUVVOaUxGTkJRV3RDTzFGQlJXeENMRWxCUVVrc1UwRkJVeXhGUVVGRk8xbEJRMklzVDBGQlR5eEpRVUZKTEVOQlFVTXNkMEpCUVhkQ0xFTkJRVU1zVVVGQlVTeEZRVUZGTEZOQlFWTXNRMEZCUXl4RFFVRkRPMU5CUXpORU8xRkJSVVFzVFVGQlRTeE5RVUZOTEVkQlFXbERMRVZCUVVVc1EwRkJRenRSUVVOb1JDeExRVUZMTEUxQlFVMHNTMEZCU3l4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVFVGQlRTeEZRVUZGTzFsQlEzWkRMRTFCUVUwc1EwRkJReXhMUVVGTExFTkJRVU1zU1VGQlNTeERRVUZETEVkQlFVY3NTVUZCU1N4RFFVRkRMSGRDUVVGM1FpeERRVUZETEZGQlFWRXNSVUZCUlN4TFFVRkxMRU5CUVVNc1NVRkJTU3hEUVVGRExFTkJRVU03VTBGRE1VVTdVVUZGUkN4UFFVRlBMRTFCUVUwc1EwRkJRenRKUVVOb1FpeERRVUZETzBsQlJVOHNkMEpCUVhkQ0xFTkJRemxDTEZGQlFXRXNSVUZEWWl4VFFVRnBRanM3VVVGRmFrSXNUVUZCVFN4WlFVRlpMRWRCUVVjc1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGRE1VSXNVMEZCUnl4RFFVRkRMRWRCUVVjc1EwRkJReXhSUVVGUkxFVkJRVVVzUTBGQlF5eEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hEUVVGRExFVkJRelZETEVsQlFVa3NRMEZCUXl4TFFVRkxMRU5CUTFnc1EwRkJRenRSUVVWR0xFMUJRVTBzUzBGQlN5eEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1RVRkJUU3hEUVVGRExFbEJRVWtzUTBGQlF5eERRVUZETEVOQlFVTXNSVUZCUlN4RlFVRkZMRU5CUVVNc1EwRkJReXhEUVVGRExFbEJRVWtzUzBGQlN5eFRRVUZUTEVOQlFVTXNRMEZCUXp0UlFVTndSU3hKUVVGSkxFTkJRVU1zUzBGQlN5eEZRVUZGTzFsQlExWXNUMEZCVHl4SlFVRkpMRU5CUVVNN1UwRkRZanRSUVVWRUxFMUJRVTBzVVVGQlVTeEhRVUZITEUxQlFVRXNUVUZCUVN4TFFVRkxMRU5CUVVNc1NVRkJTU3d3UTBGQlJ5eFpRVUZaTEVOQlFVTXNRMEZCUXl4RFFVRkRMREJEUVVGSExGbEJRVmtzUTBGQlF5eERRVUZETEVOQlFVTXNRMEZCUXp0UlFVTm9SU3hKUVVGSkxGRkJRVkVzUzBGQlN5eFRRVUZUTEVsQlFVa3NVVUZCVVN4TFFVRkxMRU5CUVVNc1EwRkJReXhGUVVGRk8xbEJRemRETEU5QlFVOHNTVUZCU1N4RFFVRkRPMU5CUTJJN1VVRkZSQ3hKUVVGSkxFdEJRVXNzUTBGQlF5eExRVUZMTEVWQlFVVTdXVUZEWml4UFFVRlBMRTFCUVVFc1MwRkJTeXhEUVVGRExFdEJRVXNzUTBGQlF5eFJRVUZSTEVOQlFVTXNiVU5CUVVrc1NVRkJTU3hEUVVGRE8xTkJRM1JETzFGQlJVUXNUMEZCVHl4SlFVRkpMRU5CUVVNN1NVRkRaQ3hEUVVGRE8wbEJSVThzVlVGQlZTeERRVUZETEVOQlFVMDdVVUZEZGtJc1QwRkJUeXhUUVVGSExFTkJRVU1zUjBGQlJ5eERRVUZETEVOQlFVTXNRMEZCUXl4RFFVRkRPMGxCUTNCQ0xFTkJRVU03U1VGRlRTeEpRVUZKTEVOQlExUXNUMEZCYVVNc1JVRkRha01zVFVGQlZ5eEZRVU5ZTEZGQlFXRXNSVUZEWWl4TFFVRmhPenRSUVVWaUxFMUJRVTBzYVVKQlFXbENMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eFJRVUZSTEVkQlFVY3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhUUVVGVExFTkJRVU03VVVGRGVrVXNUVUZCVFN4WFFVRlhMRWRCUVVjc1NVRkJRU3hUUVVGSExFVkJRVU1zU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4WFFVRlhMRU5CUVVNc1EwRkJRenRSUVVWc1JDeHZRa0ZCYjBJN1VVRkRjRUlzU1VGQlNTeFhRVUZYTEVkQlFVY3NTMEZCU3l4RFFVRkRPMUZCUTNoQ0xFbEJRVWtzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRWxCUVVrc1YwRkJWeXhIUVVGSExFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNVVUZCVVN4RlFVRkZPMWxCUTJoRkxGZEJRVmNzUjBGQlJ5eEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRkZCUVZFc1EwRkJRenRUUVVOeVF6dFJRVU5FTEVsQlFVa3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFbEJRVWtzVjBGQlZ5eEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hGUVVGRk8xbEJRMmhGTEZkQlFWY3NSMEZCUnl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGRkJRVkVzUTBGQlF6dFRRVU55UXp0UlFVVkVMR2xEUVVGcFF6dFJRVU5xUXl4SlFVRkpMR05CUVdNc1IwRkJSeXhKUVVGQkxGTkJRVWNzUlVGQlF5eFJRVUZSTEVOQlFVTXNRMEZCUXp0UlFVTnVReXhKUVVGSkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNUVUZCVFN4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zY1VKQlFYRkNMRVZCUVVVN1dVRkROMFFzVFVGQlRTeGpRVUZqTEVkQlFVY3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFZEJRVWNzVjBGQlZ5eERRVUZETzFsQlF6TkVMRTFCUVUwc1owSkJRV2RDTEVkQlFVY3NVMEZCUnl4RFFVRkRMRWRCUVVjc1EwRkRPVUlzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZCUXl4TlFVRk5MRVZCUVVVc1EwRkJReXhIUVVGSExFTkJRVU1zVjBGQlZ5eEhRVUZITEVOQlFVTXNRMEZCUXl4RFFVRkRMRVZCUTNSRExFbEJRVWtzUTBGQlF5eEpRVUZKTEVOQlExWXNRMEZCUXp0WlFVTkdMRTFCUVUwc1YwRkJWeXhIUVVGSExFbEJRVUVzVTBGQlJ5eEZRVU55UWl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFMUJRVTBzUTBGQlF5eFBRVUZQTEVOQlFVTXNRMEZCUXl4SFFVRkhMR05CUVdNc1IwRkJSeXhuUWtGQlowSXNRMEZCUXl4RFFVRkRMRVZCUTI1RkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNUVUZCVFN4RFFVRkRMRTlCUVU4c1EwRkJReXhEUVVGRExFZEJRVWNzWTBGQll5eEhRVUZITEdkQ1FVRm5RaXhEUVVGRExFTkJRVU1zUTBGRGNFVXNRMEZCUXp0WlFVTkdMRTFCUVUwc1YwRkJWeXhIUVVGSExFbEJRVUVzVTBGQlJ5eEZRVU55UWl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFMUJRVTBzUTBGQlF5eFhRVUZYTEVOQlFVTXNRMEZCUXl4SFFVRkhMR05CUVdNc1IwRkJSeXhuUWtGQlowSXNRMEZCUXl4RFFVRkRMRVZCUTNaRkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNUVUZCVFN4RFFVRkRMRmRCUVZjc1EwRkJReXhEUVVGRExFZEJRVWNzWTBGQll5eEhRVUZITEdkQ1FVRm5RaXhEUVVGRExFTkJRVU1zUTBGRGVFVXNRMEZCUXp0WlFVVkdMR05CUVdNc1IwRkJSeXhKUVVGQkxGTkJRVWNzUlVGRGJFSXNTMEZCU3l4RFFVRkRMR05CUVdNc1EwRkJReXhEUVVGRExFVkJRVVVzVjBGQlZ5eERRVUZETEVOQlFVTXNSVUZCUlN4WFFVRlhMRU5CUVVNc1EwRkJReXhEUVVGRExFVkJRM0pFTEV0QlFVc3NRMEZCUXl4alFVRmpMRU5CUVVNc1EwRkJReXhGUVVGRkxGZEJRVmNzUTBGQlF5eERRVUZETEVWQlFVVXNWMEZCVnl4RFFVRkRMRU5CUVVNc1EwRkJReXhEUVVOMFJDeERRVUZETzFOQlEwZzdVVUZGUkN4TlFVRk5MR3RDUVVGclFpeEhRVUZITEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUTJoRExGTkJRVWNzUTBGQlF5eEhRVUZITEVOQlEwd3NUVUZCVFN4RlFVTk9MRU5CUVVNc1IwRkJSeXhEUVVGRExHbENRVUZwUWl4SFFVRkhMRmRCUVZjc1EwRkJReXhEUVVOMFF5eEZRVU5FTEVsQlFVa3NRMEZCUXl4SlFVRkpMRU5CUTFZc1EwRkJRenRSUVVOR0xFMUJRVTBzYVVKQlFXbENMRWRCUVVjc1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGREwwSXNVMEZCUnl4RFFVRkRMRWRCUVVjc1EwRkJReXhqUVVGakxFVkJRVVVzUTBGQlF5eEhRVUZITEdsQ1FVRnBRaXhEUVVGRExFVkJRemxETEVsQlFVa3NRMEZCUXl4TFFVRkxMRU5CUTFnc1EwRkJRenRSUVVOR0xFMUJRVTBzV1VGQldTeEhRVUZITEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUXpGQ0xGTkJRVWNzUTBGQlF5eEhRVUZITEVOQlEwd3NhVUpCUVdsQ0xFVkJRMnBDTEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUTB3c1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGQlF5eHJRa0ZCYTBJc1JVRkJSU3hIUVVGSExFTkJRVU1zUlVGRGFFTXNTVUZCU1N4RFFVRkRMRWxCUVVrc1EwRkRWaXhEUVVOR0xFVkJRMFFzVjBGQlZ5eERRVU5hTEVOQlFVTTdVVUZEUml4TlFVRk5MR2RDUVVGblFpeEhRVUZITEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUXpsQ0xGTkJRVWNzUTBGQlF5eEhRVUZITEVOQlEwd3NhVUpCUVdsQ0xFVkJRMnBDTEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUTB3c1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGQlF5eHJRa0ZCYTBJc1JVRkJSU3hIUVVGSExFTkJRVU1zUlVGRGFFTXNTVUZCU1N4RFFVRkRMRWxCUVVrc1EwRkRWaXhEUVVOR0xFVkJRMFFzVjBGQlZ5eERRVU5hTEVOQlFVTTdVVUZGUml4UFFVRlBMRU5CUVVNc1NVRkJTU3hGUVVGRkxFTkJRVU03VVVGRFppeFBRVUZQTEVOQlFVTXNTMEZCU3l4RFFVRkRMRmRCUVZjc1JVRkJSU3hYUVVGWExFTkJRVU1zUTBGQlF6dFJRVU40UXl4UFFVRlBMRU5CUVVNc1UwRkJVeXhEUVVObUxFTkJRVU1zWTBGQll5eERRVUZETEVOQlFVTXNSMEZCUnl4TlFVRk5MRU5CUVVNc1EwRkJReXhIUVVGSExFTkJRVU1zVjBGQlZ5eEhRVUZITEVOQlFVTXNRMEZCUXl4RlFVTm9SQ3hEUVVGRExHTkJRV01zUTBGQlF5eERRVUZETEVkQlFVY3NUVUZCVFN4RFFVRkRMRU5CUVVNc1IwRkJSeXhEUVVGRExGZEJRVmNzUjBGQlJ5eERRVUZETEVOQlFVTXNRMEZEYWtRc1EwRkJRenRSUVVWR0xFMUJRVUVzVFVGQlFTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RlFVRkRMRk5CUVZNc2JVUkJRM0JDTEU5QlFVOHNSVUZEVUN4SlFVRkpMRVZCUTBvc1RVRkJUU3hGUVVOT0xHTkJRV01zUlVGRFpDeFhRVUZYTEVOQlExb3NRMEZCUXp0UlFVVkdMR2RDUVVGblFqdFJRVU5vUWl4TFFVRkxMRWxCUVVrc1EwRkJReXhIUVVGSExGbEJRVmtzUTBGQlF5eERRVUZETEVWQlFVVXNRMEZCUXl4SFFVRkhMR2RDUVVGblFpeERRVUZETEVOQlFVTXNSVUZCUlN4RFFVRkRMRVZCUVVVc1JVRkJSVHRaUVVONFJDeExRVUZMTEVsQlFVa3NRMEZCUXl4SFFVRkhMRmxCUVZrc1EwRkJReXhEUVVGRExFVkJRVVVzUTBGQlF5eEhRVUZITEdkQ1FVRm5RaXhEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVWQlFVVXNSVUZCUlR0blFrRkRlRVFzVFVGQlRTeGhRVUZoTEVkQlFVY3NTVUZCUVN4VFFVRkhMRVZCUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zUTBGQlF5eERRVUZETzJkQ1FVTm9ReXhOUVVGTkxIRkNRVUZ4UWl4SFFVRkhMRk5CUVVjc1EwRkJReXhIUVVGSExFTkJRVU1zWVVGQllTeEZRVUZGTEdsQ1FVRnBRaXhEUVVGRExFTkJRVU03WjBKQlJYaEZMREpEUVVFeVF6dG5Ra0ZETTBNc1RVRkJUU3hUUVVGVExFZEJRVWNzU1VGQlNTeERRVUZETEZWQlFWVXNRMEZCUXl4aFFVRmhMRU5CUVVNc1EwRkJRenRuUWtGRGFrUXNTVUZCU1N4RFFVRkRMRWxCUVVrc1EwRkJReXhYUVVGWExFTkJRVU1zUjBGQlJ5eERRVUZETEZOQlFWTXNRMEZCUXl4RlFVRkZPMjlDUVVOd1F5eEpRVUZKTEVOQlFVTXNWMEZCVnl4RFFVRkRMRWRCUVVjc1EwRkJReXhUUVVGVExFVkJRVVVzU1VGQlNTeERRVUZETEdGQlFXRXNRMEZEYUVRc1lVRkJZU3hGUVVOaUxHbENRVUZwUWl4RFFVTnNRaXhEUVVGRExFTkJRVU03YVVKQlEwbzdaMEpCUlVRc1RVRkJUU3hMUVVGTExFZEJRVWNzU1VGQlNTeERRVUZETEZkQlFWY3NRMEZCUXl4SFFVRkhMRU5CUVVNc1UwRkJVeXhEUVVGRExFTkJRVU03WjBKQlF6bERMRWxCUVVrc1MwRkJTeXhGUVVGRk8yOUNRVU5VTEU5QlFVOHNRMEZCUXl4VFFVRlRMRU5CUTJZc1MwRkJTeXhEUVVGRExFdEJRVXNzUlVGRFdDeHhRa0ZCY1VJc1EwRkJReXhEUVVGRExFVkJRM1pDTEhGQ1FVRnhRaXhEUVVGRExFTkJRVU1zUTBGRGVFSXNRMEZCUXp0cFFrRkRTRHRoUVVOR08xTkJRMFk3VVVGRlJDeE5RVUZCTEUxQlFVRXNTVUZCU1N4RFFVRkRMRTlCUVU4c1JVRkJReXhWUVVGVkxHMUVRVU55UWl4UFFVRlBMRVZCUTFBc1NVRkJTU3hGUVVOS0xFMUJRVTBzUlVGRFRpeGpRVUZqTEVWQlEyUXNWMEZCVnl4RFFVTmFMRU5CUVVNN1VVRkZSaXgxUWtGQmRVSTdVVUZEZGtJc1NVRkJTU3hKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEV0QlFVc3NRMEZCUXl4bFFVRmxMRVZCUVVVN1dVRkRkRU1zVFVGQlRTeFhRVUZYTEVkQlFVY3NVMEZCUnl4RFFVRkRMRWRCUVVjc1EwRkRla0lzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZEVEN4cFFrRkJhVUlzUlVGRGFrSXNVMEZCUnl4RFFVRkRMRWRCUVVjc1EwRkRUQ3hUUVVGSExFTkJRVU1zUjBGQlJ5eERRVU5NTEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUVVNc2EwSkJRV3RDTEVWQlFVVXNSMEZCUnl4RFFVRkRMRVZCUTJoRExFbEJRVWtzUTBGQlF5eEpRVUZKTEVOQlExWXNSVUZEUkN4SlFVRkJMRk5CUVVjc1JVRkJReXhEUVVGRExFTkJRVU1zUTBGRFVDeERRVU5HTEVWQlEwUXNTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhUUVVGVExFTkJRM1pDTEVOQlFVTTdXVUZEUml4TlFVRk5MR1ZCUVdVc1IwRkJSeXhUUVVGSExFTkJRVU1zUjBGQlJ5eERRVU0zUWl4VFFVRkhMRU5CUVVNc1IwRkJSeXhEUVVOTUxHbENRVUZwUWl4RlFVTnFRaXhUUVVGSExFTkJRVU1zUjBGQlJ5eERRVU5NTEZOQlFVY3NRMEZCUXl4SFFVRkhMRU5CUTB3c1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGQlF5eHJRa0ZCYTBJc1JVRkJSU3hIUVVGSExFTkJRVU1zUlVGRGFFTXNTVUZCU1N4RFFVRkRMRWxCUVVrc1EwRkRWaXhGUVVORUxFbEJRVUVzVTBGQlJ5eEZRVUZETEVOQlFVTXNRMEZCUXl4RFFVTlFMRU5CUTBZc1JVRkRSQ3hKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEZOQlFWTXNRMEZEZGtJc1EwRkJRenRaUVVWR0xFdEJRVXNzU1VGQlNTeERRVUZETEVkQlFVY3NWMEZCVnl4RFFVRkRMRU5CUVVNc1JVRkJSU3hEUVVGRExFZEJRVWNzWlVGQlpTeERRVUZETEVOQlFVTXNSVUZCUlN4RFFVRkRMRVZCUVVVc1JVRkJSVHRuUWtGRGRFUXNTVUZCU1N4RFFVRkRMRkZCUVZFc1EwRkRXQ3hQUVVGUExFVkJRMUFzU1VGQlFTeFRRVUZITEVWQlEwUXNZMEZCWXl4RFFVRkRMRU5CUVVNc1IwRkJSeXhOUVVGTkxFTkJRVU1zUTBGQlF5eEhRVUZITEVOQlFVTXNWMEZCVnl4SFFVRkhMRU5CUVVNc1EwRkJReXhGUVVNdlF5eERRVUZETEVkQlFVY3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFTkJRekZDTEVWQlEwUXNTVUZCUVN4VFFVRkhMRVZCUTBRc1kwRkJZeXhEUVVGRExFTkJRVU1zUjBGQlJ5eE5RVUZOTEVOQlFVTXNRMEZCUXl4SFFVRkhMRU5CUVVNc1YwRkJWeXhIUVVGSExFTkJRVU1zUTBGQlF5eEZRVU12UXl4RFFVRkRMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eFJRVUZSTEVOQlF6RkNMRVZCUTBRc1QwRkJUeXhEUVVGRExIZENRVUYzUWl4RlFVTm9ReXhQUVVGUExFTkJRVU1zTkVKQlFUUkNMRU5CUTNKRExFTkJRVU03WVVGRFNEdFpRVU5FTEV0QlFVc3NTVUZCU1N4RFFVRkRMRWRCUVVjc1YwRkJWeXhEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVkQlFVY3NaVUZCWlN4RFFVRkRMRU5CUVVNc1JVRkJSU3hEUVVGRExFVkJRVVVzUlVGQlJUdG5Ra0ZEZEVRc1NVRkJTU3hEUVVGRExGRkJRVkVzUTBGRFdDeFBRVUZQTEVWQlExQXNTVUZCUVN4VFFVRkhMRVZCUTBRc1EwRkJReXhIUVVGSExFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNVVUZCVVN4RlFVTjZRaXhqUVVGakxFTkJRVU1zUTBGQlF5eEhRVUZITEUxQlFVMHNRMEZCUXl4RFFVRkRMRWRCUVVjc1EwRkJReXhYUVVGWExFZEJRVWNzUTBGQlF5eERRVUZETEVOQlEyaEVMRVZCUTBRc1NVRkJRU3hUUVVGSExFVkJRMFFzUTBGQlF5eEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hGUVVONlFpeGpRVUZqTEVOQlFVTXNRMEZCUXl4SFFVRkhMRTFCUVUwc1EwRkJReXhEUVVGRExFZEJRVWNzUTBGQlF5eFhRVUZYTEVkQlFVY3NRMEZCUXl4RFFVRkRMRU5CUTJoRUxFVkJRMFFzVDBGQlR5eERRVUZETEhkQ1FVRjNRaXhGUVVOb1F5eFBRVUZQTEVOQlFVTXNORUpCUVRSQ0xFTkJRM0pETEVOQlFVTTdZVUZEU0R0VFFVTkdPMUZCUlVRc1NVRkJTU3hKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEV0QlFVc3NRMEZCUXl4blFrRkJaMElzUlVGQlJUdFpRVU4yUXl4TFFVRkxMRWxCUVVrc1EwRkJReXhIUVVGSExGbEJRVmtzUTBGQlF5eERRVUZETEVWQlFVVXNRMEZCUXl4SFFVRkhMR2RDUVVGblFpeERRVUZETEVOQlFVTXNSVUZCUlN4RFFVRkRMRVZCUVVVc1JVRkJSVHRuUWtGRGVFUXNTVUZCU1N4RFFVRkRMRkZCUVZFc1EwRkRXQ3hQUVVGUExFVkJRMUFzU1VGQlFTeFRRVUZITEVWQlEwUXNZMEZCWXl4RFFVRkRMRU5CUVVNc1IwRkJSeXhOUVVGTkxFTkJRVU1zUTBGQlF5eEhRVUZITEVOQlFVTXNWMEZCVnl4SFFVRkhMRU5CUVVNc1EwRkJReXhGUVVNdlF5eERRVUZETEVkQlFVY3NhVUpCUVdsQ0xFTkJRM1JDTEVWQlEwUXNTVUZCUVN4VFFVRkhMRVZCUTBRc1kwRkJZeXhEUVVGRExFTkJRVU1zUjBGQlJ5eE5RVUZOTEVOQlFVTXNRMEZCUXl4SFFVRkhMRU5CUVVNc1YwRkJWeXhIUVVGSExFTkJRVU1zUTBGQlF5eEZRVU12UXl4RFFVRkRMRWRCUVVjc2FVSkJRV2xDTEVOQlEzUkNMRVZCUTBRc1QwRkJUeXhEUVVGRExIbENRVUY1UWl4RlFVTnFReXhQUVVGUExFTkJRVU1zTmtKQlFUWkNMRU5CUTNSRExFTkJRVU03WVVGRFNEdFpRVU5FTEV0QlFVc3NTVUZCU1N4RFFVRkRMRWRCUVVjc1dVRkJXU3hEUVVGRExFTkJRVU1zUlVGQlJTeERRVUZETEVkQlFVY3NaMEpCUVdkQ0xFTkJRVU1zUTBGQlF5eEZRVUZGTEVOQlFVTXNSVUZCUlN4RlFVRkZPMmRDUVVONFJDeEpRVUZKTEVOQlFVTXNVVUZCVVN4RFFVTllMRTlCUVU4c1JVRkRVQ3hKUVVGQkxGTkJRVWNzUlVGRFJDeERRVUZETEVkQlFVY3NhVUpCUVdsQ0xFVkJRM0pDTEdOQlFXTXNRMEZCUXl4RFFVRkRMRWRCUVVjc1RVRkJUU3hEUVVGRExFTkJRVU1zUjBGQlJ5eERRVUZETEZkQlFWY3NSMEZCUnl4RFFVRkRMRU5CUVVNc1EwRkRhRVFzUlVGRFJDeEpRVUZCTEZOQlFVY3NSVUZEUkN4RFFVRkRMRWRCUVVjc2FVSkJRV2xDTEVWQlEzSkNMR05CUVdNc1EwRkJReXhEUVVGRExFZEJRVWNzVFVGQlRTeERRVUZETEVOQlFVTXNSMEZCUnl4RFFVRkRMRmRCUVZjc1IwRkJSeXhEUVVGRExFTkJRVU1zUTBGRGFFUXNSVUZEUkN4UFFVRlBMRU5CUVVNc2VVSkJRWGxDTEVWQlEycERMRTlCUVU4c1EwRkJReXcyUWtGQk5rSXNRMEZEZEVNc1EwRkJRenRoUVVOSU8xTkJRMFk3VVVGRlJDeEpRVUZKTEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1MwRkJTeXhEUVVGRExHVkJRV1VzUlVGQlJUdFpRVU4wUXl4UFFVRlBMRU5CUVVNc1NVRkJTU3hGUVVGRkxFTkJRVU03V1VGRFppeFBRVUZQTEVOQlFVTXNVMEZCVXl4SFFVRkhMRTlCUVU4c1EwRkJReXgzUWtGQmQwSXNRMEZCUXp0WlFVTnlSQ3hQUVVGUExFTkJRVU1zU1VGQlNTeEhRVUZITEU5QlFVOHNRMEZCUXl4elFrRkJjMElzUTBGQlF6dFpRVU01UXl4UFFVRlBMRU5CUVVNc1dVRkJXU3hIUVVGSExGRkJRVkVzUTBGQlF6dFpRVU5vUXl4UFFVRlBMRU5CUVVNc1UwRkJVeXhIUVVGSExGRkJRVkVzUTBGQlF6dFpRVVUzUWl4TFFVRkxMRWxCUVVrc1EwRkJReXhIUVVGSExGbEJRVmtzUTBGQlF5eERRVUZETEVWQlFVVXNRMEZCUXl4SFFVRkhMR2RDUVVGblFpeERRVUZETEVOQlFVTXNSVUZCUlN4RFFVRkRMRVZCUVVVc1JVRkJSVHRuUWtGRGVFUXNTMEZCU3l4SlFVRkpMRU5CUVVNc1IwRkJSeXhaUVVGWkxFTkJRVU1zUTBGQlF5eEZRVUZGTEVOQlFVTXNSMEZCUnl4blFrRkJaMElzUTBGQlF5eERRVUZETEVWQlFVVXNRMEZCUXl4RlFVRkZMRVZCUVVVN2IwSkJRM2hFTEU5QlFVOHNRMEZCUXl4UlFVRlJMRU5CUTJRc1IwRkJSeXhEUVVGRExFdEJRVXNzUTBGQlF5eEZRVUZGTEVWQlExb3NRMEZCUXl4SFFVRkhMR2xDUVVGcFFpeEhRVUZITEdsQ1FVRnBRaXhIUVVGSExFTkJRVU1zUlVGRE4wTXNRMEZCUXl4SFFVRkhMR2xDUVVGcFFpeEhRVUZITEdsQ1FVRnBRaXhIUVVGSExFTkJRVU1zUTBGRE9VTXNRMEZCUXp0cFFrRkRTRHRoUVVOR08xbEJSVVFzVDBGQlR5eERRVUZETEU5QlFVOHNSVUZCUlN4RFFVRkRPMU5CUTI1Q08xRkJSVVFzU1VGRFJTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRXRCUVVzc1EwRkJReXhWUVVGVk8xbEJRemRDTEdkQ1FVRm5RaXhEUVVGRExFbEJRVUVzVTBGQlJ5eEZRVUZETEVOQlFVTXNSVUZCUlN4RFFVRkRMRU5CUVVNc1JVRkJSU3haUVVGWkxFVkJRVVVzWjBKQlFXZENMRU5CUVVNc1JVRkRNMFE3V1VGRFFTeEpRVUZKTEVOQlFVTXNVMEZCVXl4RFFVTmFMRTlCUVU4c1JVRkRVQ3hKUVVGQkxGTkJRVWNzUlVGQlF5eERRVUZETEVWQlFVVXNRMEZCUXl4RFFVRkRMRVZCUTFRc1QwRkJUeXhEUVVGRExHMUNRVUZ0UWl4RlFVTXpRaXhQUVVGUExFTkJRVU1zZFVKQlFYVkNMRVZCUXk5Q0xFOUJRVThzUTBGQlF5eHBRa0ZCYVVJc1EwRkRNVUlzUTBGQlF6dFRRVU5JTzFGQlJVUXNUMEZCVHl4RFFVRkRMRTlCUVU4c1JVRkJSU3hEUVVGRE8wbEJRM0JDTEVOQlFVTTdTVUZGVHl4aFFVRmhMRU5CUTI1Q0xHRkJRV3RDTEVWQlEyeENMR2xDUVVGNVFqczdVVUZGZWtJc1RVRkJUU3hYUVVGWExFZEJRVWNzVVVGQlVTeERRVUZETEdGQlFXRXNRMEZCUXl4UlFVRlJMRU5CUVVNc1EwRkJRenRSUVVOeVJDeE5RVUZOTEZsQlFWa3NSMEZCUnl4WFFVRlhMRU5CUVVNc1ZVRkJWU3hEUVVGRExFbEJRVWtzUTBGQlJTeERRVUZETzFGQlJXNUVMRmRCUVZjc1EwRkJReXhMUVVGTExFZEJRVWNzYVVKQlFXbENMRU5CUVVNN1VVRkRkRU1zVjBGQlZ5eERRVUZETEUxQlFVMHNSMEZCUnl4cFFrRkJhVUlzUTBGQlF6dFJRVVYyUXl4SlFVRkpMRXRCUVVzc1IwRkJhVUk3V1VGRGVFSXNZVUZCWVR0WlFVTmlMRXRCUVVzc1JVRkJSU3hYUVVGWE8xTkJRMjVDTEVOQlFVTTdVVUZGUml4TlFVRk5MRmRCUVZjc1IwRkJSeXhUUVVGSExFTkJRVU1zUjBGQlJ5eERRVUZETEdGQlFXRXNSVUZCUlN4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGTkJRVk1zUTBGQlF5eERRVUZETzFGQlEyNUZMRTFCUVUwc1pVRkJaU3hIUVVGSExGTkJRVWNzUTBGQlF5eEhRVUZITEVOQlF6ZENMRmRCUVZjc1JVRkRXQ3hKUVVGQkxGTkJRVWNzUlVGQlF5eEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRk5CUVZNc1IwRkJSeXhEUVVGRExFTkJRVU1zUTBGRGFFTXNRMEZCUXp0UlFVTkdMRTFCUVUwc1lVRkJZU3hIUVVGSExFMUJRVUVzVFVGQlFTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRTFCUVUwc01FTkJRVVVzVDBGQlR5eHRRMEZCU1N4SlFVRkJMRk5CUVVjc1JVRkJReXhEUVVGRExFTkJRVU1zUTBGQlF6dFJRVVUzUkN4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zWjBKQlFXZENMRVZCUVVVN1dVRkRha01zVFVGQlRTeE5RVUZOTEVkQlFVY3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhuUWtGQlowSXNRMEZETVVNc1dVRkJXU3hGUVVOYUxFbEJRVWtzUlVGRFNqdG5Ra0ZEUlN4UFFVRlBMRVZCUVVVc1YwRkJWenRuUWtGRGNFSXNWMEZCVnl4RlFVRkZMR1ZCUVdVN1lVRkROMElzUlVGRFJDeGhRVUZoTEVOQlEyUXNRMEZCUXp0WlFVVkdMRWxCUVVrc1MwRkJTeXhEUVVGRExFOUJRVThzUTBGQlF5eE5RVUZOTEVOQlFVTXNSVUZCUlR0blFrRkRla0lzU1VGQlNTeERRVUZETEUxQlFVMHNRMEZCUXl4RFFVRkRMRU5CUVVNc1JVRkJSVHR2UWtGRFpDeFBRVUZQTEV0QlFVc3NRMEZCUXp0cFFrRkRaRHRoUVVOR08xTkJRMFk3VVVGRlJDeHhSRUZCY1VRN1VVRkRja1FzUzBGQlN5eE5RVUZOTEV0QlFVc3NTVUZCU1N4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFMUJRVTBzUlVGQlJUdFpRVU4yUXl4WlFVRlpMRU5CUVVNc1NVRkJTU3hGUVVGRkxFTkJRVU03V1VGRGNFSXNXVUZCV1N4RFFVRkRMRmRCUVZjc1IwRkJSeXhOUVVGQkxFdEJRVXNzUTBGQlF5eFBRVUZQTEcxRFFVRkpMRU5CUVVNc1EwRkJRenRaUVVVNVF5eE5RVUZOTEZOQlFWTXNSMEZCUnl4TlFVRkJMRXRCUVVzc1EwRkJReXhUUVVGVExHMURRVUZKTEdGQlFXRXNRMEZCUXl4TlFVRk5MRU5CUVVNN1dVRkZNVVFzUzBGQlN5eEpRVUZKTEVOQlFVTXNSMEZCUnl4WFFVRlhMRU5CUVVNc1EwRkJReXhGUVVGRkxFTkJRVU1zU1VGQlNTeGxRVUZsTEVOQlFVTXNRMEZCUXl4RlFVRkZMRU5CUVVNc1JVRkJSU3hGUVVGRk8yZENRVU4yUkN4TFFVRkxMRWxCUVVrc1EwRkJReXhIUVVGSExGZEJRVmNzUTBGQlF5eERRVUZETEVWQlFVVXNRMEZCUXl4SlFVRkpMR1ZCUVdVc1EwRkJReXhEUVVGRExFVkJRVVVzUTBGQlF5eEZRVUZGTEVWQlFVVTdiMEpCUTNaRUxFMUJRVTBzV1VGQldTeEhRVUZITEVsQlFVRXNVMEZCUnl4RlFVRkRMRU5CUVVNc1JVRkJSU3hEUVVGRExFTkJRVU1zUTBGQlF6dHZRa0ZGTDBJc1RVRkJRU3hMUVVGTExFTkJRVU1zWVVGQllTeHpSRUZEYWtJc1dVRkJXU3hGUVVOYUxFbEJRVWtzUlVGRFNpeExRVUZMTEVWQlEwd3NZVUZCWVN4RlFVTmlMRmxCUVZrc1EwRkRZaXhEUVVGRE8yOUNRVVZHTEUxQlFVMHNaMEpCUVdkQ0xFZEJRVWNzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZET1VJc1dVRkJXU3hGUVVOYUxHRkJRV0VzUTBGRFpDeERRVUZETzI5Q1FVVkdMRWxCUVVrc1owSkJRV2RDTEVOQlFVTXNRMEZCUXl4SFFVRkhMRU5CUVVNc1NVRkJTU3huUWtGQlowSXNRMEZCUXl4RFFVRkRMRWRCUVVjc1EwRkJReXhGUVVGRk8zZENRVU53UkN4VFFVRlRPM0ZDUVVOV08yOUNRVVZFTEUxQlFVMHNVVUZCVVN4SFFVRkhMRTFCUVVFc1RVRkJRU3hMUVVGTExFTkJRVU1zU1VGQlNTd3dRMEZEZEVJc1owSkJRV2RDTEVOQlFVTXNRMEZCUXl4RFFVRkRMREJEUVVOdVFpeG5Ra0ZCWjBJc1EwRkJReXhEUVVGRExFTkJRVU1zUTBGQlF6dHZRa0ZEZWtJc1NVRkJTU3hSUVVGUkxFdEJRVXNzVTBGQlV5eEpRVUZKTEZGQlFWRXNTMEZCU3l4RFFVRkRMRU5CUVVNc1JVRkJSVHQzUWtGRE4wTXNVMEZCVXp0eFFrRkRWanR2UWtGRlJDeE5RVUZOTEZOQlFWTXNSMEZCUnl4TlFVRkJMRTFCUVVFc1MwRkJTeXhEUVVGRExFdEJRVXNzTUVOQlFVY3NVVUZCVVN4RFFVRkRMREJEUVVGRkxFdEJRVXNzUTBGQlF6dHZRa0ZEYWtRc1NVRkJTU3hEUVVGRExGTkJRVk1zUlVGQlJUdDNRa0ZEWkN4VFFVRlRPM0ZDUVVOV08yOUNRVVZFTEUxQlFVMHNiMEpCUVc5Q0xFZEJRVWNzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZEYkVNc1UwRkJSeXhEUVVGRExFZEJRVWNzUTBGRFRDeFpRVUZaTEVWQlExb3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFTkJRM1JDTEVWQlEwUXNVMEZCUnl4RFFVRkRMRWRCUVVjc1EwRkJReXhoUVVGaExFVkJRVVVzYVVKQlFXbENMRU5CUVVNc1EwRkRNVU1zUTBGQlF6dHZRa0ZGUml4blFrRkJaMEk3YjBKQlEyaENMRWxCUVVrc1MwRkJTeXhEUVVGRExFbEJRVWtzUlVGQlJUdDNRa0ZEWkN4WlFVRlpMRU5CUVVNc1NVRkJTU3hGUVVGRkxFTkJRVU03ZDBKQlEzQkNMRmxCUVZrc1EwRkJReXhUUVVGVExFVkJRVVVzUTBGQlF6dDNRa0ZEZWtJc1dVRkJXU3hEUVVGRExFbEJRVWtzUTBGRFppeHZRa0ZCYjBJc1EwRkJReXhEUVVGRExFVkJRM1JDTEc5Q1FVRnZRaXhEUVVGRExFTkJRVU1zUlVGRGRFSXNTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFVkJRM0pDTEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hEUVVOMFFpeERRVUZETzNkQ1FVTkdMRmxCUVZrc1EwRkJReXhKUVVGSkxFVkJRVVVzUTBGQlF6dHhRa0ZEY2tJN2IwSkJSVVFzYVVKQlFXbENPMjlDUVVOcVFpeEpRVUZKTEhsQ1FVRTRRaXhEUVVGRE8yOUNRVU51UXl4UlFVRlJMRk5CUVZNc1JVRkJSVHQzUWtGRGFrSXNTMEZCU3l4aFFVRmhMRU5CUVVNc1QwRkJUenMwUWtGRGVFSXNlVUpCUVhsQ0xFZEJRVWNzU1VGQlFTeFRRVUZITEVWQlFVTXNiMEpCUVc5Q0xFTkJRVU1zUTBGQlF6czBRa0ZEZEVRc1RVRkJUVHQzUWtGRlVpeExRVUZMTEdGQlFXRXNRMEZCUXl4SFFVRkhPelJDUVVOd1FpeDVRa0ZCZVVJc1IwRkJSeXhKUVVGQkxGTkJRVWNzUlVGRE4wSXNRMEZEUlN4dlFrRkJiMElzUTBGQlF5eERRVUZETEVkQlFVY3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFZEJRVWNzUTBGQlF5eERRVU51UkN4SFFVRkhMRk5CUVZNc1EwRkJReXhMUVVGTExFZEJRVWNzUTBGQlF5eEZRVU4yUWl4dlFrRkJiMElzUTBGQlF5eERRVUZETEVOQlEzWkNMRU5CUVVNN05FSkJRMFlzVFVGQlRUdDNRa0ZGVWl4TFFVRkxMR0ZCUVdFc1EwRkJReXhSUVVGUk96UkNRVU42UWl4NVFrRkJlVUlzUjBGQlJ5eEpRVUZCTEZOQlFVY3NSVUZETjBJc2IwSkJRVzlDTEVOQlFVTXNRMEZCUXl4SFFVRkhMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVVVGQlVTeEhRVUZITEZOQlFWTXNRMEZCUXl4TFFVRkxMRVZCUTJoRkxHOUNRVUZ2UWl4RFFVRkRMRU5CUVVNc1EwRkRka0lzUTBGQlF6czBRa0ZEUml4TlFVRk5PM2RDUVVWU0xFdEJRVXNzWVVGQllTeERRVUZETEVsQlFVazdORUpCUTNKQ0xIbENRVUY1UWl4SFFVRkhMRWxCUVVFc1UwRkJSeXhGUVVNM1FpeHZRa0ZCYjBJc1EwRkJReXhEUVVGRExFVkJRM1JDTEVOQlEwVXNiMEpCUVc5Q0xFTkJRVU1zUTBGQlF5eEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hIUVVGSExFTkJRVU1zUTBGRGJrUXNSMEZCUnl4VFFVRlRMRU5CUVVNc1RVRkJUU3hIUVVGSExFTkJRVU1zUTBGRGVrSXNRMEZCUXpzMFFrRkRSaXhOUVVGTk8zZENRVVZTTEV0QlFVc3NZVUZCWVN4RFFVRkRMRTFCUVUwN05FSkJRM1pDTEhsQ1FVRjVRaXhIUVVGSExFbEJRVUVzVTBGQlJ5eEZRVU0zUWl4RFFVTkZMRzlDUVVGdlFpeERRVUZETEVOQlFVTXNSMEZCUnl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGRkJRVkVzUjBGQlJ5eERRVUZETEVOQlEyNUVMRWRCUVVjc1UwRkJVeXhEUVVGRExFdEJRVXNzUjBGQlJ5eERRVUZETEVWQlEzWkNMRU5CUTBVc2IwSkJRVzlDTEVOQlFVTXNRMEZCUXl4SFFVRkhMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVVVGQlVTeEhRVUZITEVOQlFVTXNRMEZEYmtRc1IwRkJSeXhUUVVGVExFTkJRVU1zVFVGQlRTeEhRVUZITEVOQlFVTXNRMEZEZWtJc1EwRkJRenMwUWtGRFJpeE5RVUZOTzNkQ1FVVlNMRXRCUVVzc1lVRkJZU3hEUVVGRExFdEJRVXM3TkVKQlEzUkNMSGxDUVVGNVFpeEhRVUZITEVsQlFVRXNVMEZCUnl4RlFVTTNRaXh2UWtGQmIwSXNRMEZCUXl4RFFVRkRMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eFJRVUZSTEVkQlFVY3NVMEZCVXl4RFFVRkRMRXRCUVVzc1JVRkRhRVVzUTBGRFJTeHZRa0ZCYjBJc1EwRkJReXhEUVVGRExFZEJRVWNzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRWRCUVVjc1EwRkJReXhEUVVOdVJDeEhRVUZITEZOQlFWTXNRMEZCUXl4TlFVRk5MRWRCUVVjc1EwRkJReXhEUVVONlFpeERRVUZET3pSQ1FVTkdMRTFCUVUwN2QwSkJSVklzUzBGQlN5eGhRVUZoTEVOQlFVTXNWVUZCVlRzMFFrRkRNMElzZVVKQlFYbENMRWRCUVVjc1NVRkJRU3hUUVVGSExFVkJRemRDTEc5Q1FVRnZRaXhEUVVGRExFTkJRVU1zUlVGRGRFSXNiMEpCUVc5Q0xFTkJRVU1zUTBGQlF5eEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hIUVVGSExGTkJRVk1zUTBGQlF5eE5RVUZOTEVOQlEyeEZMRU5CUVVNN05FSkJRMFlzVFVGQlRUdDNRa0ZGVWl4TFFVRkxMR0ZCUVdFc1EwRkJReXhOUVVGTk96UkNRVU4yUWl4NVFrRkJlVUlzUjBGQlJ5eEpRVUZCTEZOQlFVY3NSVUZETjBJc1EwRkRSU3h2UWtGQmIwSXNRMEZCUXl4RFFVRkRMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eFJRVUZSTEVkQlFVY3NRMEZCUXl4RFFVTnVSQ3hIUVVGSExGTkJRVk1zUTBGQlF5eExRVUZMTEVkQlFVY3NRMEZCUXl4RlFVTjJRaXh2UWtGQmIwSXNRMEZCUXl4RFFVRkRMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eFJRVUZSTEVkQlFVY3NVMEZCVXl4RFFVRkRMRTFCUVUwc1EwRkRiRVVzUTBGQlF6czBRa0ZEUml4TlFVRk5PM2RDUVVWU0xFdEJRVXNzWVVGQllTeERRVUZETEZkQlFWYzdORUpCUXpWQ0xIbENRVUY1UWl4SFFVRkhMRWxCUVVFc1UwRkJSeXhGUVVNM1FpeHZRa0ZCYjBJc1EwRkJReXhEUVVGRExFZEJRVWNzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRWRCUVVjc1UwRkJVeXhEUVVGRExFdEJRVXNzUlVGRGFFVXNiMEpCUVc5Q0xFTkJRVU1zUTBGQlF5eEhRVUZITEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hIUVVGSExGTkJRVk1zUTBGQlF5eE5RVUZOTEVOQlEyeEZMRU5CUVVNN05FSkJRMFlzVFVGQlRUdHhRa0ZEVkR0dlFrRkZSQ3haUVVGWkxFTkJRVU1zVTBGQlV5eERRVU53UWl4VFFVRlRMRVZCUTFRc2VVSkJRWGxDTEVOQlFVTXNRMEZCUXl4RlFVTXpRaXg1UWtGQmVVSXNRMEZCUXl4RFFVRkRMRU5CUXpWQ0xFTkJRVU03YjBKQlJVWXNTVUZCU1N4TFFVRkxMRU5CUVVNc1NVRkJTU3hGUVVGRk8zZENRVU5rTEZsQlFWa3NRMEZCUXl4UFFVRlBMRVZCUVVVc1EwRkJRenR4UWtGRGVFSTdiMEpCUlVRc1RVRkJRU3hMUVVGTExFTkJRVU1zWTBGQll5eHpSRUZEYkVJc1YwRkJWeXhGUVVOWUxGbEJRVmtzUlVGRFdpeEpRVUZKTEVWQlEwb3NTMEZCU3l4RlFVTk1MR0ZCUVdFc1JVRkRZaXhaUVVGWkxFTkJRMklzUTBGQlF6dHBRa0ZEU0R0aFFVTkdPMWxCUlVRc1dVRkJXU3hEUVVGRExFOUJRVThzUlVGQlJTeERRVUZETzFOQlEzaENPMUZCUlVRc1RVRkJRU3hOUVVGQkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVWQlFVTXNhVUpCUVdsQ0xHMUVRVU0xUWl4WFFVRlhMRVZCUTFnc1dVRkJXU3hGUVVOYUxFbEJRVWtzUlVGRFNqdFpRVU5GTEU5QlFVOHNSVUZCUlN4WFFVRlhPMWxCUTNCQ0xGZEJRVmNzUlVGQlJTeGxRVUZsTzFOQlF6ZENMRVZCUTBRc1lVRkJZU3hEUVVOa0xFTkJRVU03VVVGRlJpeFBRVUZQTEV0QlFVc3NRMEZCUXp0SlFVTm1MRU5CUVVNN1NVRkZUeXhSUVVGUkxFTkJRMlFzVDBGQmFVTXNSVUZEYWtNc1MwRkJWU3hGUVVOV0xFZEJRVkVzUlVGRFVpeE5RVUZqTEVWQlEyUXNVMEZCYVVJN1VVRkZha0lzVDBGQlR5eERRVUZETEVsQlFVa3NSVUZCUlN4RFFVRkRPMUZCUldZc1QwRkJUeXhEUVVGRExGTkJRVk1zUjBGQlJ5eFRRVUZUTEVOQlFVTTdVVUZET1VJc1QwRkJUeXhEUVVGRExGZEJRVmNzUjBGQlJ5eE5RVUZOTEVOQlFVTTdVVUZGTjBJc1QwRkJUeXhEUVVGRExGTkJRVk1zUlVGQlJTeERRVUZETzFGQlEzQkNMRTlCUVU4c1EwRkJReXhOUVVGTkxFTkJRVU1zUzBGQlN5eERRVUZETEVOQlFVTXNSVUZCUlN4TFFVRkxMRU5CUVVNc1EwRkJReXhEUVVGRExFTkJRVU03VVVGRGFrTXNUMEZCVHl4RFFVRkRMRTFCUVUwc1EwRkJReXhIUVVGSExFTkJRVU1zUTBGQlF5eEZRVUZGTEVkQlFVY3NRMEZCUXl4RFFVRkRMRU5CUVVNc1EwRkJRenRSUVVNM1FpeFBRVUZQTEVOQlFVTXNUVUZCVFN4RlFVRkZMRU5CUVVNN1VVRkZha0lzVDBGQlR5eERRVUZETEU5QlFVOHNSVUZCUlN4RFFVRkRPMGxCUTNCQ0xFTkJRVU03U1VGRlR5eFRRVUZUTEVOQlEyWXNUMEZCYVVNc1JVRkRha01zVVVGQllTeEZRVU5pTEUxQlFXTXNSVUZEWkN4VFFVRnBRaXhGUVVOcVFpeEpRVUZaTzFGQlJWb3NUMEZCVHl4RFFVRkRMRWxCUVVrc1JVRkJSU3hEUVVGRE8xRkJSV1lzVDBGQlR5eERRVUZETEZOQlFWTXNSMEZCUnl4VFFVRlRMRU5CUVVNN1VVRkZPVUlzVFVGQlRTeFJRVUZSTEVkQlFVY3NTVUZCU1N4RFFVRkRMRWxCUVVrc1EwRkJReXhKUVVGSkxFZEJRVWNzUTBGQlF5eERRVUZETEVOQlFVTTdVVUZEY2tNc1QwRkJUeXhEUVVGRExGZEJRVmNzUjBGQlJ5eE5RVUZOTEVOQlFVTTdVVUZETjBJc1QwRkJUeXhEUVVGRExGTkJRVk1zUlVGQlJTeERRVUZETzFGQlEzQkNMRTlCUVU4c1EwRkJReXhOUVVGTkxFTkJRVU1zVVVGQlVTeERRVUZETEVOQlFVTXNSMEZCUnl4UlFVRlJMRVZCUVVVc1VVRkJVU3hEUVVGRExFTkJRVU1zUjBGQlJ5eFJRVUZSTEVOQlFVTXNRMEZCUXp0UlFVTTNSQ3hQUVVGUExFTkJRVU1zVFVGQlRTeERRVUZETEZGQlFWRXNRMEZCUXl4RFFVRkRMRWRCUVVjc1VVRkJVU3hGUVVGRkxGRkJRVkVzUTBGQlF5eERRVUZETEVkQlFVY3NVVUZCVVN4RFFVRkRMRU5CUVVNN1VVRkROMFFzVDBGQlR5eERRVUZETEUxQlFVMHNRMEZCUXl4UlFVRlJMRU5CUVVNc1EwRkJReXhIUVVGSExGRkJRVkVzUlVGQlJTeFJRVUZSTEVOQlFVTXNRMEZCUXl4SFFVRkhMRkZCUVZFc1EwRkJReXhEUVVGRE8xRkJRemRFTEU5QlFVOHNRMEZCUXl4TlFVRk5MRU5CUVVNc1VVRkJVU3hEUVVGRExFTkJRVU1zUjBGQlJ5eFJRVUZSTEVWQlFVVXNVVUZCVVN4RFFVRkRMRU5CUVVNc1IwRkJSeXhSUVVGUkxFTkJRVU1zUTBGQlF6dFJRVU0zUkN4UFFVRlBMRU5CUVVNc1RVRkJUU3hGUVVGRkxFTkJRVU03VVVGRmFrSXNUMEZCVHl4RFFVRkRMRTlCUVU4c1JVRkJSU3hEUVVGRE8wbEJRM0JDTEVOQlFVTTdPMEZCTjI1Q1NDd3dRa0U0YmtKRE8wRkJOMjVDUXl4UFFVRlBPMEZCUlZBc2QwUkJRWGRFTzBGQlEzaEVMSGxEUVVGNVF6dEJRVVZxUWl4MVFrRkJaU3hIUVVGdFFqdEpRVU40UkN4eFFrRkJjVUlzUlVGQlJTeEpRVUZKTzBsQlF6TkNMRkZCUVZFc1JVRkJSU3hGUVVGRk8wbEJRMW9zVFVGQlRTeEZRVUZGTzFGQlEwNDdXVUZEUlN4SlFVRkpMRVZCUVVVc1UwRkJVenRUUVVOb1FqdExRVU5HTzBsQlEwUXNVMEZCVXl4RlFVRkZMRU5CUVVNN1NVRkRXaXhYUVVGWExFVkJRVVVzUTBGQlF6dEpRVU5rTEd0Q1FVRnJRaXhGUVVGRkxFVkJRVVU3UTBGRGRrSXNRMEZCUXp0QlFVVnpRaXd5UWtGQmJVSXNSMEZCUnl4TlFVRk5MRU5CUVVNN1FVRkROMElzSzBKQlFYVkNMRWRCUVVjc1EwRkJReXhEUVVGRE8wRkJRelZDTEhsQ1FVRnBRaXhIUVVGSExFVkJRVVVzUTBGQlF6dEJRVVYyUWl4cFEwRkJlVUlzUjBGQlJ5eFJRVUZSTEVOQlFVTTdRVUZEY2tNc2NVTkJRVFpDTEVkQlFVY3NRMEZCUXl4RFFVRkRPMEZCUld4RExHZERRVUYzUWl4SFFVRkhMRTlCUVU4c1EwRkJRenRCUVVOdVF5dzRRa0ZCYzBJc1IwRkJSeXhuUWtGQlowSXNRMEZCUXp0QlFVVXhReXhuUTBGQmQwSXNSMEZCUnl4UlFVRlJMRU5CUVVNN1FVRkRjRU1zYjBOQlFUUkNMRWRCUVVjc1EwRkJReXhEUVVGRE8wRkJhMjFDTTBRN096czdPMGRCUzBjN1FVRkRTU3hMUVVGTExGVkJRVlVzT0VKQlFUaENMRU5CUTJ4RUxFOUJTMFVzUlVGRFJpeEpRVXRETzBsQlJVUXNSVUZCUlR0QlFVTktMRU5CUVVNN1FVRm1SQ3gzUlVGbFF5SjkiLCIvLyBUaGUgbW9kdWxlIGNhY2hlXG52YXIgX193ZWJwYWNrX21vZHVsZV9jYWNoZV9fID0ge307XG5cbi8vIFRoZSByZXF1aXJlIGZ1bmN0aW9uXG5mdW5jdGlvbiBfX3dlYnBhY2tfcmVxdWlyZV9fKG1vZHVsZUlkKSB7XG5cdC8vIENoZWNrIGlmIG1vZHVsZSBpcyBpbiBjYWNoZVxuXHR2YXIgY2FjaGVkTW9kdWxlID0gX193ZWJwYWNrX21vZHVsZV9jYWNoZV9fW21vZHVsZUlkXTtcblx0aWYgKGNhY2hlZE1vZHVsZSAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0cmV0dXJuIGNhY2hlZE1vZHVsZS5leHBvcnRzO1xuXHR9XG5cdC8vIENyZWF0ZSBhIG5ldyBtb2R1bGUgKGFuZCBwdXQgaXQgaW50byB0aGUgY2FjaGUpXG5cdHZhciBtb2R1bGUgPSBfX3dlYnBhY2tfbW9kdWxlX2NhY2hlX19bbW9kdWxlSWRdID0ge1xuXHRcdC8vIG5vIG1vZHVsZS5pZCBuZWVkZWRcblx0XHQvLyBubyBtb2R1bGUubG9hZGVkIG5lZWRlZFxuXHRcdGV4cG9ydHM6IHt9XG5cdH07XG5cblx0Ly8gRXhlY3V0ZSB0aGUgbW9kdWxlIGZ1bmN0aW9uXG5cdF9fd2VicGFja19tb2R1bGVzX19bbW9kdWxlSWRdLmNhbGwobW9kdWxlLmV4cG9ydHMsIG1vZHVsZSwgbW9kdWxlLmV4cG9ydHMsIF9fd2VicGFja19yZXF1aXJlX18pO1xuXG5cdC8vIFJldHVybiB0aGUgZXhwb3J0cyBvZiB0aGUgbW9kdWxlXG5cdHJldHVybiBtb2R1bGUuZXhwb3J0cztcbn1cblxuIiwiIiwiLy8gc3RhcnR1cFxuLy8gTG9hZCBlbnRyeSBtb2R1bGUgYW5kIHJldHVybiBleHBvcnRzXG4vLyBUaGlzIGVudHJ5IG1vZHVsZSBpcyByZWZlcmVuY2VkIGJ5IG90aGVyIG1vZHVsZXMgc28gaXQgY2FuJ3QgYmUgaW5saW5lZFxudmFyIF9fd2VicGFja19leHBvcnRzX18gPSBfX3dlYnBhY2tfcmVxdWlyZV9fKFwiLi9pbmRleC50c1wiKTtcbiIsIiJdLCJuYW1lcyI6W10sInNvdXJjZVJvb3QiOiIifQ==