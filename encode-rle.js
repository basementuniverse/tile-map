#!/usr/bin/env node

const encode = require('fast-rle/encode');

const data = JSON.parse(process.argv[2]);
const flattenedData = data.reduce((acc, val) => acc.concat(val), []);

console.log(JSON.stringify(encode(flattenedData)));
