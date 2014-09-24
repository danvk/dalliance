(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);throw new Error("Cannot find module '"+o+"'")}var f=n[o]={exports:{}};t[o][0].call(f.exports,function(e){var n=t[o][1][e];return s(n?n:e)},f,f.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bam.js: indexed binary alignments
//

"use strict";

if (typeof(require) !== 'undefined') {
    var spans = require('./spans');
    var Range = spans.Range;
    var union = spans.union;
    var intersection = spans.intersection;

    var bin = require('./bin');
    var readInt = bin.readInt;
    var readShort = bin.readShort;
    var readByte = bin.readByte;
    var readInt64 = bin.readInt64;
    var readFloat = bin.readFloat;

    var lh3utils = require('./lh3utils');
    var readVob = lh3utils.readVob;
    var unbgzf = lh3utils.unbgzf;
    var reg2bins = lh3utils.reg2bins;
    var Chunk = lh3utils.Chunk;
}


var BAM_MAGIC = 0x14d4142;
var BAI_MAGIC = 0x1494142;

var BamFlags = {
    MULTIPLE_SEGMENTS:       0x1,
    ALL_SEGMENTS_ALIGN:      0x2,
    SEGMENT_UNMAPPED:        0x4,
    NEXT_SEGMENT_UNMAPPED:   0x8,
    REVERSE_COMPLEMENT:      0x10,
    NEXT_REVERSE_COMPLEMENT: 0x20,
    FIRST_SEGMENT:           0x40,
    LAST_SEGMENT:            0x80,
    SECONDARY_ALIGNMENT:     0x100,
    QC_FAIL:                 0x200,
    DUPLICATE:               0x400,
    SUPPLEMENTARY:           0x800
};

function BamFile() {
}

function makeBam(data, bai, callback) {
    var bam = new BamFile();
    bam.data = data;
    bam.bai = bai;

    bam.bai.fetch(function(header) {   // Do we really need to fetch the whole thing? :-(
        if (!header) {
            return callback(null, "Couldn't access BAI");
        }

        var uncba = new Uint8Array(header);
        var baiMagic = readInt(uncba, 0);
        if (baiMagic != BAI_MAGIC) {
            return callback(null, 'Not a BAI file, magic=0x' + baiMagic.toString(16));
        }

        var nref = readInt(uncba, 4);

        bam.indices = [];

        var p = 8;
        var minBlockIndex = 1000000000;
        for (var ref = 0; ref < nref; ++ref) {
            var blockStart = p;
            var nbin = readInt(uncba, p); p += 4;
            for (var b = 0; b < nbin; ++b) {
                var bin = readInt(uncba, p);
                var nchnk = readInt(uncba, p+4);
                p += 8 + (nchnk * 16);
            }
            var nintv = readInt(uncba, p); p += 4;
            
            var q = p;
            for (var i = 0; i < nintv; ++i) {
                var v = readVob(uncba, q); q += 8;
                if (v) {
                    var bi = v.block;
                    if (v.offset > 0)
                        bi += 65536;

                    if (bi < minBlockIndex)
                        minBlockIndex = bi;
                    break;
                }
            }
            p += (nintv * 8);


            if (nbin > 0) {
                bam.indices[ref] = new Uint8Array(header, blockStart, p - blockStart);
            }                     
        }

        bam.data.slice(0, minBlockIndex).fetch(function(r) {
            if (!r) {
                return callback(null, "Couldn't access BAM");
            }
            
            var unc = unbgzf(r, r.byteLength);
            var uncba = new Uint8Array(unc);

            var magic = readInt(uncba, 0);
            if (magic != BAM_MAGIC) {
                return callback(null, "Not a BAM file, magic=0x" + magic.toString(16));
            }
            var headLen = readInt(uncba, 4);
            var header = '';
            for (var i = 0; i < headLen; ++i) {
                header += String.fromCharCode(uncba[i + 8]);
            }

            var nRef = readInt(uncba, headLen + 8);
            var p = headLen + 12;

            bam.chrToIndex = {};
            bam.indexToChr = [];
            for (var i = 0; i < nRef; ++i) {
                var lName = readInt(uncba, p);
                var name = '';
                for (var j = 0; j < lName-1; ++j) {
                    name += String.fromCharCode(uncba[p + 4 + j]);
                }
                var lRef = readInt(uncba, p + lName + 4);
                bam.chrToIndex[name] = i;
                if (name.indexOf('chr') == 0) {
                    bam.chrToIndex[name.substring(3)] = i;
                } else {
                    bam.chrToIndex['chr' + name] = i;
                }
                bam.indexToChr.push(name);

                p = p + 8 + lName;
            }

            if (bam.indices) {
                return callback(bam);
            }
        });
    });
}



BamFile.prototype.blocksForRange = function(refId, min, max) {
    var index = this.indices[refId];
    if (!index) {
        return [];
    }

    var intBinsL = reg2bins(min, max);
    var intBins = [];
    for (var i = 0; i < intBinsL.length; ++i) {
        intBins[intBinsL[i]] = true;
    }
    var leafChunks = [], otherChunks = [];

    var nbin = readInt(index, 0);
    var p = 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(index, p);
        var nchnk = readInt(index, p+4);
//        dlog('bin=' + bin + '; nchnk=' + nchnk);
        p += 8;
        if (intBins[bin]) {
            for (var c = 0; c < nchnk; ++c) {
                var cs = readVob(index, p);
                var ce = readVob(index, p + 8);
                (bin < 4681 ? otherChunks : leafChunks).push(new Chunk(cs, ce));
                p += 16;
            }
        } else {
            p +=  (nchnk * 16);
        }
    }
//    dlog('leafChunks = ' + miniJSONify(leafChunks));
//    dlog('otherChunks = ' + miniJSONify(otherChunks));

    var nintv = readInt(index, p);
    var lowest = null;
    var minLin = Math.min(min>>14, nintv - 1), maxLin = Math.min(max>>14, nintv - 1);
    for (var i = minLin; i <= maxLin; ++i) {
        var lb =  readVob(index, p + 4 + (i * 8));
        if (!lb) {
            continue;
        }
        if (!lowest || lb.block < lowest.block || lb.offset < lowest.offset) {
            lowest = lb;
        }
    }
    // dlog('Lowest LB = ' + lowest);
    
    var prunedOtherChunks = [];
    if (lowest != null) {
        for (var i = 0; i < otherChunks.length; ++i) {
            var chnk = otherChunks[i];
            if (chnk.maxv.block >= lowest.block && chnk.maxv.offset >= lowest.offset) {
                prunedOtherChunks.push(chnk);
            }
        }
    }
    // dlog('prunedOtherChunks = ' + miniJSONify(prunedOtherChunks));
    otherChunks = prunedOtherChunks;

    var intChunks = [];
    for (var i = 0; i < otherChunks.length; ++i) {
        intChunks.push(otherChunks[i]);
    }
    for (var i = 0; i < leafChunks.length; ++i) {
        intChunks.push(leafChunks[i]);
    }

    intChunks.sort(function(c0, c1) {
        var dif = c0.minv.block - c1.minv.block;
        if (dif != 0) {
            return dif;
        } else {
            return c0.minv.offset - c1.minv.offset;
        }
    });
    var mergedChunks = [];
    if (intChunks.length > 0) {
        var cur = intChunks[0];
        for (var i = 1; i < intChunks.length; ++i) {
            var nc = intChunks[i];
            if (nc.minv.block == cur.maxv.block /* && nc.minv.offset == cur.maxv.offset */) { // no point splitting mid-block
                cur = new Chunk(cur.minv, nc.maxv);
            } else {
                mergedChunks.push(cur);
                cur = nc;
            }
        }
        mergedChunks.push(cur);
    }
    // dlog('mergedChunks = ' + miniJSONify(mergedChunks));

    return mergedChunks;
}

BamFile.prototype.fetch = function(chr, min, max, callback, opts) {
    var thisB = this;
    opts = opts || {};

    var chrId = this.chrToIndex[chr];
    var chunks;
    if (chrId === undefined) {
        chunks = [];
    } else {
        chunks = this.blocksForRange(chrId, min, max);
        if (!chunks) {
            callback(null, 'Error in index fetch');
        }
    }
    
    var records = [];
    var index = 0;
    var data;

    function tramp() {
        if (index >= chunks.length) {
            return callback(records);
        } else if (!data) {
            // dlog('fetching ' + index);
            var c = chunks[index];
            var fetchMin = c.minv.block;
            var fetchMax = c.maxv.block + (1<<16); // *sigh*
            thisB.data.slice(fetchMin, fetchMax - fetchMin).fetch(function(r) {
                data = unbgzf(r, c.maxv.block - c.minv.block + 1);
                return tramp();
            });
        } else {
            var ba = new Uint8Array(data);
            thisB.readBamRecords(ba, chunks[index].minv.offset, records, min, max, chrId, opts);
            data = null;
            ++index;
            return tramp();
        }
    }
    tramp();
}

var SEQRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
var CIGAR_DECODER = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];

function BamRecord() {
}

BamFile.prototype.readBamRecords = function(ba, offset, sink, min, max, chrId, opts) {
    while (true) {
        var blockSize = readInt(ba, offset);
        var blockEnd = offset + blockSize + 4;
        if (blockEnd >= ba.length) {
            return sink;
        }

        var record = new BamRecord();

        var refID = readInt(ba, offset + 4);
        var pos = readInt(ba, offset + 8);
        
        var bmn = readInt(ba, offset + 12);
        var bin = (bmn & 0xffff0000) >> 16;
        var mq = (bmn & 0xff00) >> 8;
        var nl = bmn & 0xff;

        var flag_nc = readInt(ba, offset + 16);
        var flag = (flag_nc & 0xffff0000) >> 16;
        var nc = flag_nc & 0xffff;
    
        var lseq = readInt(ba, offset + 20);
        
        var nextRef  = readInt(ba, offset + 24);
        var nextPos = readInt(ba, offset + 28);
        
        var tlen = readInt(ba, offset + 32);
    
        record.segment = this.indexToChr[refID];
        record.flag = flag;
        record.pos = pos;
        record.mq = mq;
        if (opts.light)
            record.seqLength = lseq;

        if (!opts.light) {
            if (nextRef >= 0) {
                record.nextSegment = this.indexToChr[nextRef];
                record.nextPos = nextPos;
            }

            var readName = '';
            for (var j = 0; j < nl-1; ++j) {
                readName += String.fromCharCode(ba[offset + 36 + j]);
            }
            record.readName = readName;
        
            var p = offset + 36 + nl;

            var cigar = '';
            for (var c = 0; c < nc; ++c) {
                var cigop = readInt(ba, p);
                cigar = cigar + (cigop>>4) + CIGAR_DECODER[cigop & 0xf];
                p += 4;
            }
            record.cigar = cigar;
        
            var seq = '';
            var seqBytes = (lseq + 1) >> 1;
            for (var j = 0; j < seqBytes; ++j) {
                var sb = ba[p + j];
                seq += SEQRET_DECODER[(sb & 0xf0) >> 4];
                seq += SEQRET_DECODER[(sb & 0x0f)];
            }
            p += seqBytes;
            record.seq = seq;

            var qseq = '';
            for (var j = 0; j < lseq; ++j) {
                qseq += String.fromCharCode(ba[p + j] + 33);
            }
            p += lseq;
            record.quals = qseq;

            while (p < blockEnd) {
                var tag = String.fromCharCode(ba[p], ba[p + 1]);
                var type = String.fromCharCode(ba[p + 2]);
                var value;

                if (type == 'A') {
                    value = String.fromCharCode(ba[p + 3]);
                    p += 4;
                } else if (type == 'i' || type == 'I') {
                    value = readInt(ba, p + 3);
                    p += 7;
                } else if (type == 'c' || type == 'C') {
                    value = ba[p + 3];
                    p += 4;
                } else if (type == 's' || type == 'S') {
                    value = readShort(ba, p + 3);
                    p += 5;
                } else if (type == 'f') {
                    value = readFloat(ba, p + 3);
                    p += 7;
                } else if (type == 'Z' || type == 'H') {
                    p += 3;
                    value = '';
                    for (;;) {
                        var cc = ba[p++];
                        if (cc == 0) {
                            break;
                        } else {
                            value += String.fromCharCode(cc);
                        }
                    }
                } else if (type == 'B') {
                    var atype = String.fromCharCode(ba[p + 3]);
                    var alen = readInt(ba, p + 4);
                    var elen;
                    var reader;
                    if (atype == 'i' || atype == 'I' || atype == 'f') {
                        elen = 4;
                        if (atype == 'f')
                            reader = readFloat;
                        else
                            reader = readInt;
                    } else if (atype == 's' || atype == 'S') {
                        elen = 2;
                        reader = readShort;
                    } else if (atype == 'c' || atype == 'C') {
                        elen = 1;
                        reader = readByte;
                    } else {
                        throw 'Unknown array type ' + atype;
                    }

                    p += 8;
                    value = [];
                    for (var i = 0; i < alen; ++i) {
                        value.push(reader(ba, p));
                        p += elen;
                    }
                } else {
                    throw 'Unknown type '+ type;
                }
                record[tag] = value;
            }
        }

        if (!min || record.pos <= max && record.pos + lseq >= min) {
            if (chrId === undefined || refID == chrId) {
                sink.push(record);
            }
        }
        offset = blockEnd;
    }

    // Exits via top of loop.
};

if (typeof(module) !== 'undefined') {
    module.exports = {
        makeBam: makeBam,
        BAM_MAGIC: BAM_MAGIC,
        BAI_MAGIC: BAI_MAGIC,
        BamFlags: BamFlags
    };
}
},{"./bin":3,"./lh3utils":7,"./spans":9}],2:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// bigwig.js: indexed binary WIG (and BED) files
//

"use strict";


if (typeof(require) !== 'undefined') {
    var spans = require('./spans');
    var Range = spans.Range;
    var union = spans.union;
    var intersection = spans.intersection;

    var das = require('./das');
    var DASFeature = das.DASFeature;
    var DASGroup = das.DASGroup;

    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;

    var bin = require('./bin');
    var readInt = bin.readInt;

    var jszlib = require('jszlib');
    var jszlib_inflate_buffer = jszlib.inflateBuffer;
    var arrayCopy = jszlib.arrayCopy;
}

var BIG_WIG_MAGIC = 0x888FFC26;
var BIG_WIG_MAGIC_BE = 0x26FC8F88;
var BIG_BED_MAGIC = 0x8789F2EB;
var BIG_BED_MAGIC_BE = 0xEBF28987;


var BIG_WIG_TYPE_GRAPH = 1;
var BIG_WIG_TYPE_VSTEP = 2;
var BIG_WIG_TYPE_FSTEP = 3;
  
var M1 = 256;
var M2 = 256*256;
var M3 = 256*256*256;
var M4 = 256*256*256*256;

var BED_COLOR_REGEXP = new RegExp("^[0-9]+,[0-9]+,[0-9]+");

function bwg_readOffset(ba, o) {
    var offset = ba[o] + ba[o+1]*M1 + ba[o+2]*M2 + ba[o+3]*M3 + ba[o+4]*M4;
    return offset;
}

function BigWig() {
}

BigWig.prototype.readChromTree = function(callback) {
    var thisB = this;
    this.chromsToIDs = {};
    this.idsToChroms = {};
    this.maxID = 0;

    var udo = this.unzoomedDataOffset;
    var eb = (udo - this.chromTreeOffset) & 3;
    udo = udo + 4 - eb;

    this.data.slice(this.chromTreeOffset, udo - this.chromTreeOffset).fetch(function(bpt) {
        var ba = new Uint8Array(bpt);
        var sa = new Int16Array(bpt);
        var la = new Int32Array(bpt);
        var bptMagic = la[0];
        var blockSize = la[1];
        var keySize = la[2];
        var valSize = la[3];
        var itemCount = bwg_readOffset(ba, 16);
        var rootNodeOffset = 32;

        var bptReadNode = function(offset) {
            var nodeType = ba[offset];
            var cnt = sa[(offset/2) + 1];
            offset += 4;
            for (var n = 0; n < cnt; ++n) {
                if (nodeType == 0) {
                    offset += keySize;
                    var childOffset = bwg_readOffset(ba, offset);
                    offset += 8;
                    childOffset -= thisB.chromTreeOffset;
                    bptReadNode(childOffset);
                } else {
                    var key = '';
                    for (var ki = 0; ki < keySize; ++ki) {
                        var charCode = ba[offset++];
                        if (charCode != 0) {
                            key += String.fromCharCode(charCode);
                        }
                    }
                    var chromId = (ba[offset+3]<<24) | (ba[offset+2]<<16) | (ba[offset+1]<<8) | (ba[offset+0]);
                    var chromSize = (ba[offset + 7]<<24) | (ba[offset+6]<<16) | (ba[offset+5]<<8) | (ba[offset+4]);
                    offset += 8;

                    thisB.chromsToIDs[key] = chromId;
                    if (key.indexOf('chr') == 0) {
                        thisB.chromsToIDs[key.substr(3)] = chromId;
                    }
                    thisB.idsToChroms[chromId] = key;
                    thisB.maxID = Math.max(thisB.maxID, chromId);
                }
            }
        };
        bptReadNode(rootNodeOffset);

        callback(thisB);
    });
}

function BigWigView(bwg, cirTreeOffset, cirTreeLength, isSummary) {
    this.bwg = bwg;
    this.cirTreeOffset = cirTreeOffset;
    this.cirTreeLength = cirTreeLength;
    this.isSummary = isSummary;
}



BigWigView.prototype.readWigData = function(chrName, min, max, callback) {
    var chr = this.bwg.chromsToIDs[chrName];
    if (chr === undefined) {
        // Not an error because some .bwgs won't have data for all chromosomes.
        return callback([]);
    } else {
        this.readWigDataById(chr, min, max, callback);
    }
}

BigWigView.prototype.readWigDataById = function(chr, min, max, callback) {
    var thisB = this;
    if (!this.cirHeader) {
        this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function(result) {
            thisB.cirHeader = result;
            var la = new Int32Array(thisB.cirHeader);
            thisB.cirBlockSize = la[1];
            thisB.readWigDataById(chr, min, max, callback);
        });
        return;
    }

    var blocksToFetch = [];
    var outstanding = 0;

    var beforeBWG = Date.now();

    var filter = function(chromId, fmin, fmax, toks) {
        return ((chr < 0 || chromId == chr) && fmin <= max && fmax >= min);
    }

    var cirFobRecur = function(offset, level) {
        if (thisB.bwg.instrument)
            console.log('level=' + level + '; offset=' + offset + '; time=' + (Date.now()|0));

        outstanding += offset.length;

        if (offset.length == 1 && offset[0] - thisB.cirTreeOffset == 48 && thisB.cachedCirRoot) {
            cirFobRecur2(thisB.cachedCirRoot, 0, level);
            --outstanding;
            if (outstanding == 0) {
                thisB.fetchFeatures(filter, blocksToFetch, callback);
            }
            return;
        }

        var maxCirBlockSpan = 4 +  (thisB.cirBlockSize * 32);   // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], offset[i] + maxCirBlockSpan);
            spans = spans ? union(spans, blockSpan) : blockSpan;
        }
        
        var fetchRanges = spans.ranges();
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            cirFobStartFetch(offset, fr, level);
        }
    }

    var cirFobStartFetch = function(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        thisB.bwg.data.slice(fr.min(), fr.max() - fr.min()).fetch(function(resultBuffer) {
            for (var i = 0; i < offset.length; ++i) {
                if (fr.contains(offset[i])) {
                    cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);

                    if (offset[i] - thisB.cirTreeOffset == 48 && offset[i] - fr.min() == 0)
                        thisB.cachedCirRoot = resultBuffer;

                    --outstanding;
                    if (outstanding == 0) {
                        thisB.fetchFeatures(filter, blocksToFetch, callback);
                    }
                }
            }
        });
    }

    var cirFobRecur2 = function(cirBlockData, offset, level) {
        var ba = new Uint8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData);
        var la = new Int32Array(cirBlockData);

        var isLeaf = ba[offset];
        var cnt = sa[offset/2 + 1];
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                var blockSize = bwg_readOffset(ba, offset+24);
                if (((chr < 0 || startChrom < chr) || (startChrom == chr && startBase <= max)) &&
                    ((chr < 0 || endChrom   > chr) || (endChrom == chr && endBase >= min)))
                {
                    blocksToFetch.push({offset: blockOffset, size: blockSize});
                }
                offset += 32;
            }
        } else {
            var recurOffsets = [];
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                if ((chr < 0 || startChrom < chr || (startChrom == chr && startBase <= max)) &&
                    (chr < 0 || endChrom   > chr || (endChrom == chr && endBase >= min)))
                {
                    recurOffsets.push(blockOffset);
                }
                offset += 24;
            }
            if (recurOffsets.length > 0) {
                cirFobRecur(recurOffsets, level + 1);
            }
        }
    };

    cirFobRecur([thisB.cirTreeOffset + 48], 1);
}


BigWigView.prototype.fetchFeatures = function(filter, blocksToFetch, callback) {
    var thisB = this;

    blocksToFetch.sort(function(b0, b1) {
        return (b0.offset|0) - (b1.offset|0);
    });

    if (blocksToFetch.length == 0) {
        callback([]);
    } else {
        var features = [];
        var createFeature = function(chr, fmin, fmax, opts) {
            if (!opts) {
                opts = {};
            }
        
            var f = new DASFeature();
            f._chromId = chr;
            f.segment = thisB.bwg.idsToChroms[chr];
            f.min = fmin;
            f.max = fmax;
            f.type = 'bigwig';
            
            for (var k in opts) {
                f[k] = opts[k];
            }
            
            features.push(f);
        };

        var tramp = function() {
            if (blocksToFetch.length == 0) {
                var afterBWG = Date.now();
                // dlog('BWG fetch took ' + (afterBWG - beforeBWG) + 'ms');
                callback(features);
                return;  // just in case...
            } else {
                var block = blocksToFetch[0];
                if (block.data) {
                    thisB.parseFeatures(block.data, createFeature, filter);
                    blocksToFetch.splice(0, 1);
                    tramp();
                } else {
                    var fetchStart = block.offset;
                    var fetchSize = block.size;
                    var bi = 1;
                    while (bi < blocksToFetch.length && blocksToFetch[bi].offset == (fetchStart + fetchSize)) {
                        fetchSize += blocksToFetch[bi].size;
                        ++bi;
                    }

                    thisB.bwg.data.slice(fetchStart, fetchSize).fetch(function(result) {
                        var offset = 0;
                        var bi = 0;
                        while (offset < fetchSize) {
                            var fb = blocksToFetch[bi];
                        
                            var data;
                            if (thisB.bwg.uncompressBufSize > 0) {
                                data = jszlib_inflate_buffer(result, offset + 2, fb.size - 2);
                            } else {
                                var tmp = new Uint8Array(fb.size);    // FIXME is this really the best we can do?
                                arrayCopy(new Uint8Array(result, offset, fb.size), 0, tmp, 0, fb.size);
                                data = tmp.buffer;
                            }
                            fb.data = data;
                            
                            offset += fb.size;
                            ++bi;
                        }
                        tramp();
                    });
                }
            }
        }
        tramp();
    }
}

BigWigView.prototype.parseFeatures = function(data, createFeature, filter) {
    var ba = new Uint8Array(data);

    if (this.isSummary) {
        var sa = new Int16Array(data);
        var la = new Int32Array(data);
        var fa = new Float32Array(data);

        var itemCount = data.byteLength/32;
        for (var i = 0; i < itemCount; ++i) {
            var chromId =   la[(i*8)];
            var start =     la[(i*8)+1];
            var end =       la[(i*8)+2];
            var validCnt =  la[(i*8)+3];
            var minVal    = fa[(i*8)+4];
            var maxVal    = fa[(i*8)+5];
            var sumData   = fa[(i*8)+6];
            var sumSqData = fa[(i*8)+7];
            
            if (filter(chromId, start + 1, end)) {
                var summaryOpts = {type: 'bigwig', score: sumData/validCnt, maxScore: maxVal};
                if (this.bwg.type == 'bigbed') {
                    summaryOpts.type = 'density';
                }
                createFeature(chromId, start + 1, end, summaryOpts);
            }
        }
    } else if (this.bwg.type == 'bigwig') {
        var sa = new Int16Array(data);
        var la = new Int32Array(data);
        var fa = new Float32Array(data);

        var chromId = la[0];
        var blockStart = la[1];
        var blockEnd = la[2];
        var itemStep = la[3];
        var itemSpan = la[4];
        var blockType = ba[20];
        var itemCount = sa[11];
        
        if (blockType == BIG_WIG_TYPE_FSTEP) {
            for (var i = 0; i < itemCount; ++i) {
                var score = fa[i + 6];
                var fmin = blockStart + (i*itemStep) + 1, fmax = blockStart + (i*itemStep) + itemSpan;
                if (filter(chromId, fmin, fmax))
                    createFeature(chromId, fmin, fmax, {score: score});
            }
        } else if (blockType == BIG_WIG_TYPE_VSTEP) {
            for (var i = 0; i < itemCount; ++i) {
                var start = la[(i*2) + 6] + 1;
                var end = start + itemSpan - 1;
                var score = fa[(i*2) + 7];
                if (filter(chromId, start, end))
                    createFeature(chromId, start, end, {score: score});
            }
        } else if (blockType == BIG_WIG_TYPE_GRAPH) {
            for (var i = 0; i < itemCount; ++i) {
                var start = la[(i*3) + 6] + 1;
                var end   = la[(i*3) + 7];
                var score = fa[(i*3) + 8];
                if (start > end) {
                    start = end;
                }
                if (filter(chromId, start, end))
                    createFeature(chromId, start, end, {score: score});
            }
        } else {
            console.log('Currently not handling bwgType=' + blockType);
        }
    } else if (this.bwg.type == 'bigbed') {
        var offset = 0;
        var dfc = this.bwg.definedFieldCount;
        var schema = this.bwg.schema;

        while (offset < ba.length) {
            var chromId = (ba[offset+3]<<24) | (ba[offset+2]<<16) | (ba[offset+1]<<8) | (ba[offset+0]);
            var start = (ba[offset+7]<<24) | (ba[offset+6]<<16) | (ba[offset+5]<<8) | (ba[offset+4]);
            var end = (ba[offset+11]<<24) | (ba[offset+10]<<16) | (ba[offset+9]<<8) | (ba[offset+8]);
            offset += 12;
            var rest = '';
            while (true) {
                var ch = ba[offset++];
                if (ch != 0) {
                    rest += String.fromCharCode(ch);
                } else {
                    break;
                }
            }

            var featureOpts = {};
            
            var bedColumns;
            if (rest.length > 0) {
                bedColumns = rest.split('\t');
            } else {
                bedColumns = [];
            }
            if (bedColumns.length > 0 && dfc > 3) {
                featureOpts.label = bedColumns[0];
            }
            if (bedColumns.length > 1 && dfc > 4) {
                var score = parseInt(bedColumns[1]);
                if (!isNaN(score))
                    featureOpts.score = score;
            }
            if (bedColumns.length > 2 && dfc > 5) {
                featureOpts.orientation = bedColumns[2];
            }
            if (bedColumns.length > 5 && dfc > 8) {
                var color = bedColumns[5];
                if (BED_COLOR_REGEXP.test(color)) {
                    featureOpts.itemRgb = 'rgb(' + color + ')';
                }
            }

            if (bedColumns.length > dfc-3 && schema) {
                for (var col = dfc - 3; col < bedColumns.length; ++col) {
                    featureOpts[schema.fields[col+3].name] = bedColumns[col];
                }
            }

            if (filter(chromId, start + 1, end, bedColumns)) {
                if (dfc < 12) {
                    createFeature(chromId, start + 1, end, featureOpts);
                } else {
                    var thickStart = bedColumns[3]|0;
                    var thickEnd   = bedColumns[4]|0;
                    var blockCount = bedColumns[6]|0;
                    var blockSizes = bedColumns[7].split(',');
                    var blockStarts = bedColumns[8].split(',');

                    if (featureOpts.exonFrames) {
                        var exonFrames = featureOpts.exonFrames.split(',');
                        featureOpts.exonFrames = undefined;
                    }
                    
                    featureOpts.type = 'transcript'
                    var grp = new DASGroup();
                    for (var k in featureOpts) {
                        grp[k] = featureOpts[k];
                    }
                    grp.id = bedColumns[0];
                    grp.segment = this.bwg.idsToChroms[chromId];
                    grp.min = start + 1;
                    grp.max = end;
                    grp.notes = [];
                    featureOpts.groups = [grp];

                    // Moving towards using bigGenePred model, but will
                    // still support old Dalliance-style BED12+gene-name for the
                    // foreseeable future.
                    if (bedColumns.length > 9) {
                        var geneId = featureOpts.geneName || bedColumns[9];
                        var geneName = geneId;
                        if (bedColumns.length > 10) {
                            geneName = bedColumns[10];
                        }
                        if (featureOpts.geneName2)
                            geneName = featureOpts.geneName2;

                        var gg = shallowCopy(grp);
                        gg.id = geneId;
                        gg.label = geneName;
                        gg.type = 'gene';
                        featureOpts.groups.push(gg);
                    }

                    var spanList = [];
                    for (var b = 0; b < blockCount; ++b) {
                        var bmin = (blockStarts[b]|0) + start;
                        var bmax = bmin + (blockSizes[b]|0);
                        var span = new Range(bmin, bmax);
                        spanList.push(span);
                    }
                    var spans = union(spanList);
                    
                    var tsList = spans.ranges();
                    for (var s = 0; s < tsList.length; ++s) {
                        var ts = tsList[s];
                        createFeature(chromId, ts.min() + 1, ts.max(), featureOpts);
                    }

                    if (thickEnd > thickStart) {
                        var codingRegion = (featureOpts.orientation == '+') ?
                            new Range(thickStart, thickEnd + 3) :
                            new Range(thickStart - 3, thickEnd);
                            // +/- 3 to account for stop codon

                        var tl = intersection(spans, codingRegion);
                        if (tl) {
                            featureOpts.type = 'translation';
                            var tlList = tl.ranges();
                            var readingFrame = 0;

                            var tlOffset = 0;
                            while (tlList[0].min() > tsList[tlOffset].max())
                                tlOffset++;

                            for (var s = 0; s < tlList.length; ++s) {
                                // Record reading frame for every exon
                                var index = s;
                                if (featureOpts.orientation == '-')
                                    index = tlList.length - s - 1;
                                var ts = tlList[index];
                                featureOpts.readframe = readingFrame;
                                if (exonFrames) {
                                    var brf = parseInt(exonFrames[index + tlOffset]);
                                    if (typeof(brf) === 'number' && brf >= 0 && brf <= 2) {
                                        featureOpts.readframe = brf;
                                        featureOpts.readframeExplicit = true;
                                    }
                                }
                                var length = ts.max() - ts.min();
                                readingFrame = (readingFrame + length) % 3;
                                createFeature(chromId, ts.min() + 1, ts.max(), featureOpts);
                            }
                        }
                    }
                }
            }
        }
    } else {
        throw Error("Don't know what to do with " + this.bwg.type);
    }
}

//
// nasty cut/paste, should roll back in!
//

BigWigView.prototype.getFirstAdjacent = function(chrName, pos, dir, callback) {
    var chr = this.bwg.chromsToIDs[chrName];
    if (chr === undefined) {
        // Not an error because some .bwgs won't have data for all chromosomes.
        return callback([]);
    } else {
        this.getFirstAdjacentById(chr, pos, dir, callback);
    }
}

BigWigView.prototype.getFirstAdjacentById = function(chr, pos, dir, callback) {
    var thisB = this;
    if (!this.cirHeader) {
        this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function(result) {
            thisB.cirHeader = result;
            var la = new Int32Array(thisB.cirHeader);
            thisB.cirBlockSize = la[1];
            thisB.getFirstAdjacentById(chr, pos, dir, callback);
        });
        return;
    }

    var blockToFetch = null;
    var bestBlockChr = -1;
    var bestBlockOffset = -1;

    var outstanding = 0;

    var beforeBWG = Date.now();

    var cirFobRecur = function(offset, level) {
        outstanding += offset.length;

        var maxCirBlockSpan = 4 +  (thisB.cirBlockSize * 32);   // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], offset[i] + maxCirBlockSpan);
            spans = spans ? union(spans, blockSpan) : blockSpan;
        }
        
        var fetchRanges = spans.ranges();
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            cirFobStartFetch(offset, fr, level);
        }
    }

    var cirFobStartFetch = function(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        thisB.bwg.data.slice(fr.min(), fr.max() - fr.min()).fetch(function(resultBuffer) {
            for (var i = 0; i < offset.length; ++i) {
                if (fr.contains(offset[i])) {
                    cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);
                    --outstanding;
                    if (outstanding == 0) {
                        if (!blockToFetch) {
                            if (dir > 0 && (chr != 0 || pos > 0)) {
                                return thisB.getFirstAdjacentById(0, 0, dir, callback);
                            } else if (dir < 0 && (chr != thisB.bwg.maxID || pos < 1000000000)) {
                                return thisB.getFirstAdjacentById(thisB.bwg.maxID, 1000000000, dir, callback);
                            }
                            return callback([]);
                        }

                        thisB.fetchFeatures(function(chrx, fmin, fmax, toks) {
                            return (dir < 0 && (chrx < chr || fmax < pos)) || (dir > 0 && (chrx > chr || fmin > pos));
                        }, [blockToFetch], function(features) {
                            var bestFeature = null;
                            var bestChr = -1;
                            var bestPos = -1;
                            for (var fi = 0; fi < features.length; ++fi) {
                                var f = features[fi];
                                var chrx = f._chromId, fmin = f.min, fmax = f.max;
                                if (bestFeature == null || ((dir < 0) && (chrx > bestChr || fmax > bestPos)) || ((dir > 0) && (chrx < bestChr || fmin < bestPos))) {
                                    bestFeature = f;
                                    bestPos = (dir < 0) ? fmax : fmin;
                                    bestChr = chrx;
                                }
                            }

                            if (bestFeature != null) 
                                return callback([bestFeature]);
                            else
                                return callback([]);
                        });
                    }
                }
            }
        });
    }

    var cirFobRecur2 = function(cirBlockData, offset, level) {
        var ba = new Uint8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData);
        var la = new Int32Array(cirBlockData);

        var isLeaf = ba[offset];
        var cnt = sa[offset/2 + 1];
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                var blockSize = bwg_readOffset(ba, offset+24);
                if ((dir < 0 && ((startChrom < chr || (startChrom == chr && startBase <= pos)))) ||
                    (dir > 0 && ((endChrom > chr || (endChrom == chr && endBase >= pos)))))
                {
                    // console.log('Got an interesting block: startBase=' + startChrom + ':' + startBase + '; endBase=' + endChrom + ':' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                    if (/_random/.exec(thisB.bwg.idsToChroms[startChrom])) {
                        // dlog('skipping random: ' + thisB.bwg.idsToChroms[startChrom]);
                    } else if (blockToFetch == null || ((dir < 0) && (endChrom > bestBlockChr || (endChrom == bestBlockChr && endBase > bestBlockOffset)) ||
                                                 (dir > 0) && (startChrom < bestBlockChr || (startChrom == bestBlockChr && startBase < bestBlockOffset))))
                    {
                        //                        dlog('best is: startBase=' + startChrom + ':' + startBase + '; endBase=' + endChrom + ':' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                        blockToFetch = {offset: blockOffset, size: blockSize};
                        bestBlockOffset = (dir < 0) ? endBase : startBase;
                        bestBlockChr = (dir < 0) ? endChrom : startChrom;
                    }
                }
                offset += 32;
            }
        } else {
            var bestRecur = -1;
            var bestPos = -1;
            var bestChr = -1;
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = (la[lo + 4]<<32) | (la[lo + 5]);
                if ((dir < 0 && ((startChrom < chr || (startChrom == chr && startBase <= pos)) &&
                                 (endChrom   >= chr))) ||
                     (dir > 0 && ((endChrom > chr || (endChrom == chr && endBase >= pos)) &&
                                  (startChrom <= chr))))
                {
                    if (bestRecur < 0 || endBase > bestPos) {
                        bestRecur = blockOffset;
                        bestPos = (dir < 0) ? endBase : startBase;
                        bestChr = (dir < 0) ? endChrom : startChrom;
                    }
                }
                offset += 24;
            }
            if (bestRecur >= 0) {
                cirFobRecur([bestRecur], level + 1);
            }
        }
    };
    

    cirFobRecur([thisB.cirTreeOffset + 48], 1);
}

BigWig.prototype.readWigData = function(chrName, min, max, callback) {
    this.getUnzoomedView().readWigData(chrName, min, max, callback);
}

BigWig.prototype.getUnzoomedView = function() {
    if (!this.unzoomedView) {
        var cirLen = 4000;
        var nzl = this.zoomLevels[0];
        if (nzl) {
            cirLen = this.zoomLevels[0].dataOffset - this.unzoomedIndexOffset;
        }
        this.unzoomedView = new BigWigView(this, this.unzoomedIndexOffset, cirLen, false);
    }
    return this.unzoomedView;
}

BigWig.prototype.getZoomedView = function(z) {
    var zh = this.zoomLevels[z];
    if (!zh.view) {
        zh.view = new BigWigView(this, zh.indexOffset, /* this.zoomLevels[z + 1].dataOffset - zh.indexOffset */ 4000, true);
    }
    return zh.view;
}

function makeBwg(data, callback, name) {
    var bwg = new BigWig();
    bwg.data = data;
    bwg.name = name;
    bwg.data.slice(0, 512).salted().fetch(function(result) {
        if (!result) {
            return callback(null, "Couldn't fetch file");
        }

        var header = result;
        var ba = new Uint8Array(header);
        var sa = new Int16Array(header);
        var la = new Int32Array(header);
        var magic = ba[0] + (M1 * ba[1]) + (M2 * ba[2]) + (M3 * ba[3]);
        if (magic == BIG_WIG_MAGIC) {
            bwg.type = 'bigwig';
        } else if (magic == BIG_BED_MAGIC) {
            bwg.type = 'bigbed';
        } else if (magic == BIG_WIG_MAGIC_BE || magic == BIG_BED_MAGIC_BE) {
            callback(null, "Currently don't support big-endian BBI files");
        } else {
            callback(null, "Not a supported format, magic=0x" + magic.toString(16));
        }

        bwg.version = sa[2];             // 4
        bwg.numZoomLevels = sa[3];       // 6
        bwg.chromTreeOffset = bwg_readOffset(ba, 8);
        bwg.unzoomedDataOffset = bwg_readOffset(ba, 16);
        bwg.unzoomedIndexOffset = bwg_readOffset(ba, 24);
        bwg.fieldCount = sa[16];         // 32
        bwg.definedFieldCount = sa[17];  // 34
        bwg.asOffset = bwg_readOffset(ba, 36);
        bwg.totalSummaryOffset = bwg_readOffset(ba, 44);
        bwg.uncompressBufSize = la[13];  // 52
        bwg.extHeaderOffset = bwg_readOffset(ba, 56);

        bwg.zoomLevels = [];
        for (var zl = 0; zl < bwg.numZoomLevels; ++zl) {
            var zlReduction = la[zl*6 + 16]
            var zlData = bwg_readOffset(ba, zl*24 + 72);
            var zlIndex = bwg_readOffset(ba, zl*24 + 80);
            bwg.zoomLevels.push({reduction: zlReduction, dataOffset: zlData, indexOffset: zlIndex});
        }

        bwg.readChromTree(function() {
            bwg.getAutoSQL(function(as) {
                bwg.schema = as;
                return callback(bwg);
            });
        });
    });
}


BigWig.prototype._tsFetch = function(zoom, chr, min, max, callback) {
    var bwg = this;
    if (zoom >= this.zoomLevels.length - 1) {
        if (!this.topLevelReductionCache) {
            this.getZoomedView(this.zoomLevels.length - 1).readWigDataById(-1, 0, 300000000, function(feats) {
                bwg.topLevelReductionCache = feats;
                return bwg._tsFetch(zoom, chr, min, max, callback);
            });
        } else {
            var f = [];
            var c = this.topLevelReductionCache;
            for (var fi = 0; fi < c.length; ++fi) {
                if (c[fi]._chromId == chr) {
                    f.push(c[fi]);
                }
            }
            return callback(f);
        }
    } else {
        var view;
        if (zoom < 0) {
            view = this.getUnzoomedView();
        } else {
            view = this.getZoomedView(zoom);
        }
        return view.readWigDataById(chr, min, max, callback);
    }
}

BigWig.prototype.thresholdSearch = function(chrName, referencePoint, dir, threshold, callback) {
    dir = (dir<0) ? -1 : 1;
    var bwg = this;
    var initialChr = this.chromsToIDs[chrName];
    var candidates = [{chrOrd: 0, chr: initialChr, zoom: bwg.zoomLevels.length - 4, min: 0, max: 300000000, fromRef: true}]
    for (var i = 1; i <= this.maxID + 1; ++i) {
        var chrId = (initialChr + (dir*i)) % (this.maxID + 1);
        if (chrId < 0) 
            chrId += (this.maxID + 1);
        candidates.push({chrOrd: i, chr: chrId, zoom: bwg.zoomLevels.length - 1, min: 0, max: 300000000})
    }
       
    function fbThresholdSearchRecur() {
    	if (candidates.length == 0) {
    	    return callback(null);
    	}
    	candidates.sort(function(c1, c2) {
    	    var d = c1.zoom - c2.zoom;
    	    if (d != 0)
    		    return d;

            d = c1.chrOrd - c2.chrOrd;
            if (d != 0)
                return d;
    	    else
    		    return c1.min - c2.min * dir;
    	});

	    var candidate = candidates.splice(0, 1)[0];
        bwg._tsFetch(candidate.zoom, candidate.chr, candidate.min, candidate.max, function(feats) {
            var rp = dir > 0 ? 0 : 300000000;
            if (candidate.fromRef)
                rp = referencePoint;
            
            for (var fi = 0; fi < feats.length; ++fi) {
    	        var f = feats[fi];
                var score;
                if (f.maxScore != undefined)
                    score = f.maxScore;
                else
                    score = f.score;

                if (dir > 0) {
    	            if (score > threshold) {
        		        if (candidate.zoom < 0) {
        		            if (f.min > rp)
                                return callback(f);
        		        } else if (f.max > rp) {
        		            candidates.push({chr: candidate.chr, chrOrd: candidate.chrOrd, zoom: candidate.zoom - 2, min: f.min, max: f.max, fromRef: candidate.fromRef});
        		        }
                    }
                } else {
                    if (score > threshold) {
            		    if (candidate.zoom < 0) {
                	        if (f.max < rp)
                			    return callback(f);
                        } else if (f.min < rp) {
                            candidates.push({chr: candidate.chr, chrOrd: candidate.chrOrd, zoom: candidate.zoom - 2, min: f.min, max: f.max, fromRef: candidate.fromRef});
                        }
    	            }
                }
    	    }
            fbThresholdSearchRecur();
        });
    }
    
    fbThresholdSearchRecur();
}

BigWig.prototype.getAutoSQL = function(callback) {
    var thisB = this;
    if (!this.asOffset)
        return callback(null);


    this.data.slice(this.asOffset, 2048).fetch(function(result) {
        var ba = new Uint8Array(result);
        var s = '';
        for (var i = 0; i < ba.length; ++i) {
            if (ba[i] == 0)
                break;
            s += String.fromCharCode(ba[i]);
        }
        
        /* 
         * Quick'n'dirty attempt to parse autoSql format.
         * See: http://www.linuxjournal.com/files/linuxjournal.com/linuxjournal/articles/059/5949/5949l2.html
         */

        var header_re = /(\w+)\s+(\w+)\s+("([^"]+)")?\s+\(\s*/;
        var field_re = /([\w\[\]]+)\s+(\w+)\s*;\s*("([^"]+)")?\s*/g;

        var headerMatch = header_re.exec(s);
        if (headerMatch) {
            var as = {
                declType: headerMatch[1],
                name: headerMatch[2],
                comment: headerMatch[4],

                fields: []
            };

            s = s.substring(headerMatch[0]);
            for (var m = field_re.exec(s); m != null; m = field_re.exec(s)) {
                as.fields.push({type: m[1],
                             name: m[2],
                             comment: m[4]});
            }

            return callback(as);
        }
    });
}

BigWig.prototype.getExtraIndices = function(callback) {
    var thisB = this;
    if (this.version < 4 || this.extHeaderOffset == 0 || this.type != 'bigbed') {
        return callback(null);
    } else {
        this.data.slice(this.extHeaderOffset, 64).fetch(function(result) {
            if (!result) {
                return callback(null, "Couldn't fetch extension header");
            }

            var ba = new Uint8Array(result);
            var sa = new Int16Array(result);
            var la = new Int32Array(result);
            
            var extHeaderSize = sa[0];
            var extraIndexCount = sa[1];
            var extraIndexListOffset = bwg_readOffset(ba, 4);

            if (extraIndexCount == 0) {
                return callback(null);
            }

            // FIXME 20byte records only make sense for single-field indices.
            // Right now, these seem to be the only things around, but the format
            // is actually more general.
            thisB.data.slice(extraIndexListOffset, extraIndexCount * 20).fetch(function(eil) {
                if (!eil) {
                    return callback(null, "Couldn't fetch index info");
                }

                var ba = new Uint8Array(eil);
                var sa = new Int16Array(eil);
                var la = new Int32Array(eil);

                var indices = [];
                for (var ii = 0; ii < extraIndexCount; ++ii) {
                    var eiType = sa[ii*10];
                    var eiFieldCount = sa[ii*10 + 1];
                    var eiOffset = bwg_readOffset(ba, ii*20 + 4);
                    var eiField = sa[ii*10 + 8]
                    var index = new BBIExtraIndex(thisB, eiType, eiFieldCount, eiOffset, eiField);
                    indices.push(index);
                }
                callback(indices);
            });
        });
    }
}

function BBIExtraIndex(bbi, type, fieldCount, offset, field) {
    this.bbi = bbi;
    this.type = type;
    this.fieldCount = fieldCount;
    this.offset = offset;
    this.field = field;
}

BBIExtraIndex.prototype.lookup = function(name, callback) {
    var thisB = this;

    this.bbi.data.slice(this.offset, 32).fetch(function(bpt) {
        var ba = new Uint8Array(bpt);
        var sa = new Int16Array(bpt);
        var la = new Int32Array(bpt);
        var bptMagic = la[0];
        var blockSize = la[1];
        var keySize = la[2];
        var valSize = la[3];
        var itemCount = bwg_readOffset(ba, 16);
        var rootNodeOffset = 32;

        function bptReadNode(nodeOffset) {
            thisB.bbi.data.slice(nodeOffset, 4 + (blockSize * (keySize + valSize))).fetch(function(node) {
                var ba = new Uint8Array(node);
                var sa = new Uint16Array(node);
                var la = new Uint32Array(node);

                var nodeType = ba[0];
                var cnt = sa[1];

                var offset = 4;
                if (nodeType == 0) {
                    var lastChildOffset = null;
                    for (var n = 0; n < cnt; ++n) {
                        var key = '';
                        for (var ki = 0; ki < keySize; ++ki) {
                            var charCode = ba[offset++];
                            if (charCode != 0) {
                                key += String.fromCharCode(charCode);
                            }
                        }

                        var childOffset = bwg_readOffset(ba, offset);
                        offset += 8;
                        
                        if (name.localeCompare(key) < 0 && lastChildOffset) {
                            bptReadNode(lastChildOffset);
                            return;
                        }
                        lastChildOffset = childOffset;
                    }
                    bptReadNode(lastChildOffset);
                } else {
                    for (var n = 0; n < cnt; ++n) {
                        var key = '';
                        for (var ki = 0; ki < keySize; ++ki) {
                            var charCode = ba[offset++];
                            if (charCode != 0) {
                                key += String.fromCharCode(charCode);
                            }
                        }
                        
                        // Specific for EI case.
                        if (key == name) {
                            var start = bwg_readOffset(ba, offset);
                            var length = readInt(ba, offset + 8);

                            return thisB.bbi.getUnzoomedView().fetchFeatures(
                                function(chr, min, max, toks) {
                                    if (toks && toks.length > thisB.field - 3)
                                        return toks[thisB.field - 3] == name;
                                }, 
                                [{offset: start, size: length}], 
                                callback);
                        }
                        offset += valSize;
                    }
                    return callback([]);
                }
            });
        }

        bptReadNode(thisB.offset + rootNodeOffset);
    });
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        makeBwg: makeBwg,
        BIG_BED_MAGIC: BIG_BED_MAGIC,
        BIG_WIG_MAGIC: BIG_WIG_MAGIC
    }
}

},{"./bin":3,"./das":5,"./spans":9,"./utils":10,"jszlib":11}],3:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bin.js general binary data support
//

"use strict";

if (typeof(require) !== 'undefined') {
    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;

    var sha1 = require('./sha1');
    var b64_sha1 = sha1.b64_sha1;
}

function BlobFetchable(b) {
    this.blob = b;
}

BlobFetchable.prototype.slice = function(start, length) {
    var b;

    if (this.blob.slice) {
        if (length) {
            b = this.blob.slice(start, start + length);
        } else {
            b = this.blob.slice(start);
        }
    } else {
        if (length) {
            b = this.blob.webkitSlice(start, start + length);
        } else {
            b = this.blob.webkitSlice(start);
        }
    }
    return new BlobFetchable(b);
}

BlobFetchable.prototype.salted = function() {return this;}

if (typeof(FileReader) !== 'undefined') {
    // console.log('defining async BlobFetchable.fetch');

    BlobFetchable.prototype.fetch = function(callback) {
        var reader = new FileReader();
        reader.onloadend = function(ev) {
            callback(bstringToBuffer(reader.result));
        };
        reader.readAsBinaryString(this.blob);
    }

} else {
    // if (console && console.log)
    //    console.log('defining sync BlobFetchable.fetch');

    BlobFetchable.prototype.fetch = function(callback) {
        var reader = new FileReaderSync();
        try {
            var res = reader.readAsArrayBuffer(this.blob);
            callback(res);
        } catch (e) {
            callback(null, e);
        }
    }
}

function URLFetchable(url, start, end, opts) {
    if (!opts) {
        if (typeof start === 'object') {
            opts = start;
            start = undefined;
        } else {
            opts = {};
        }
    }

    this.url = url;
    this.start = start || 0;
    if (end) {
        this.end = end;
    }
    this.opts = opts;
}

URLFetchable.prototype.slice = function(s, l) {
    if (s < 0) {
        throw 'Bad slice ' + s;
    }

    var ns = this.start, ne = this.end;
    if (ns && s) {
        ns = ns + s;
    } else {
        ns = s || ns;
    }
    if (l && ns) {
        ne = ns + l - 1;
    } else {
        ne = ne || l - 1;
    }
    return new URLFetchable(this.url, ns, ne, this.opts);
}

var seed=0;
var isSafari = navigator.userAgent.indexOf('Safari') >= 0 && navigator.userAgent.indexOf('Chrome') < 0 ;

URLFetchable.prototype.fetchAsText = function(callback) {
    var req = new XMLHttpRequest();
    var length;
    var url = this.url;
    if (isSafari || this.opts.salt) {
        url = url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++seed));
    }
    req.open('GET', url, true);

    if (this.end) {
        if (this.end - this.start > 100000000) {
            throw 'Monster fetch!';
        }
        req.setRequestHeader('Range', 'bytes=' + this.start + '-' + this.end);
        length = this.end - this.start + 1;
    }

    req.onreadystatechange = function() {
        if (req.readyState == 4) {
            if (req.status == 200 || req.status == 206) {
                return callback(req.responseText);
            } else {
                return callback(null);
            }
        }
    };
    if (this.opts.credentials) {
        req.withCredentials = true;
    }
    req.send('');
}

URLFetchable.prototype.salted = function() {
    var o = shallowCopy(this.opts);
    o.salt = true;
    return new URLFetchable(this.url, this.start, this.end, o);
}

URLFetchable.prototype.fetch = function(callback, attempt, truncatedLength) {
    var thisB = this;

    attempt = attempt || 1;
    if (attempt > 3) {
        return callback(null);
    }

    var req = new XMLHttpRequest();
    var length;
    var url = this.url;
    if (isSafari || this.opts.salt) {
        url = url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++seed));
    }
    req.open('GET', url, true);
    req.overrideMimeType('text/plain; charset=x-user-defined');
    if (this.end) {
        if (this.end - this.start > 100000000) {
            throw 'Monster fetch!';
        }
        req.setRequestHeader('Range', 'bytes=' + this.start + '-' + this.end);
        length = this.end - this.start + 1;
    }
    req.responseType = 'arraybuffer';
    req.onreadystatechange = function() {
        if (req.readyState == 4) {
            if (req.status == 200 || req.status == 206) {
                if (req.response) {
                    var bl = req.response.byteLength;
                    if (length && length != bl && (!truncatedLength || bl != truncatedLength)) {
                        return thisB.fetch(callback, attempt + 1, bl);
                    } else {
                        return callback(req.response);
                    }
                } else if (req.mozResponseArrayBuffer) {
                    return callback(req.mozResponseArrayBuffer);
                } else {
                    var r = req.responseText;
                    if (length && length != r.length && (!truncatedLength || r.length != truncatedLength)) {
                        return thisB.fetch(callback, attempt + 1, r.length);
                    } else {
                        return callback(bstringToBuffer(req.responseText));
                    }
                }
            } else {
                return thisB.fetch(callback, attempt + 1);
            }
        }
    };
    if (this.opts.credentials) {
        req.withCredentials = true;
    }
    req.send('');
}

function bstringToBuffer(result) {
    if (!result) {
        return null;
    }

    var ba = new Uint8Array(result.length);
    for (var i = 0; i < ba.length; ++i) {
        ba[i] = result.charCodeAt(i);
    }
    return ba.buffer;
}

// Read from Uint8Array

(function(global) {
    var convertBuffer = new ArrayBuffer(8);
    var ba = new Uint8Array(convertBuffer);
    var fa = new Float32Array(convertBuffer);


    global.readFloat = function(buf, offset) {
        ba[0] = buf[offset];
        ba[1] = buf[offset+1];
        ba[2] = buf[offset+2];
        ba[3] = buf[offset+3];
        return fa[0];
    };
 }(this));

function readInt64(ba, offset) {
    return (ba[offset + 7] << 24) | (ba[offset + 6] << 16) | (ba[offset + 5] << 8) | (ba[offset + 4]);
}

function readInt(ba, offset) {
    return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
}

function readShort(ba, offset) {
    return (ba[offset + 1] << 8) | (ba[offset]);
}

function readByte(ba, offset) {
    return ba[offset];
}

function readIntBE(ba, offset) {
    return (ba[offset] << 24) | (ba[offset + 1] << 16) | (ba[offset + 2] << 8) | (ba[offset + 3]);
}

// Exports if we are being used as a module

if (typeof(module) !== 'undefined') {
    module.exports = {
        BlobFetchable: BlobFetchable,
        URLFetchable: URLFetchable,

        readInt: readInt,
        readIntBE: readIntBE,
        readInt64: readInt64,
        readShort: readShort,
        readByte: readByte,
        readFloat: this.readFloat
    }
}

},{"./sha1":8,"./utils":10}],4:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// color.js
//

"use strict";

function DColour(red, green, blue, name) {
    this.red = red|0;
    this.green = green|0;
    this.blue = blue|0;
    if (name) {
        this.name = name;
    }
}

DColour.prototype.toSvgString = function() {
    if (!this.name) {
        this.name = "rgb(" + this.red + "," + this.green + "," + this.blue + ")";
    }

    return this.name;
}

function hex2(x) {
    var y = '00' + x.toString(16);
    return y.substring(y.length - 2);
}

DColour.prototype.toHexString = function() {
    return '#' + hex2(this.red) + hex2(this.green) + hex2(this.blue);
}

var palette = {
    red: new DColour(255, 0, 0, 'red'),
    green: new DColour(0, 255, 0, 'green'),
    blue: new DColour(0, 0, 255, 'blue'),
    yellow: new DColour(255, 255, 0, 'yellow'),
    white: new DColour(255, 255, 255, 'white'),
    black: new DColour(0, 0, 0, 'black'),
    gray: new DColour(180, 180, 180, 'gray'),
    grey: new DColour(180, 180, 180, 'grey'),
    lightskyblue: new DColour(135, 206, 250, 'lightskyblue'),
    lightsalmon: new DColour(255, 160, 122, 'lightsalmon'),
    hotpink: new DColour(255, 105, 180, 'hotpink')
};

var COLOR_RE = new RegExp('^#([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})$');
var CSS_COLOR_RE = /rgb\(([0-9]+),([0-9]+),([0-9]+)\)/

function dasColourForName(name) {
    var c = palette[name];
    if (!c) {
        var match = COLOR_RE.exec(name);
        if (match) {
            c = new DColour(('0x' + match[1])|0, ('0x' + match[2])|0, ('0x' + match[3])|0, name);
            palette[name] = c;
        } else {
    	    match = CSS_COLOR_RE.exec(name);
    	    if (match) {
        		c = new DColour(match[1]|0, match[2]|0, match[3]|0, name);
        		palette[name] = c;
	       } else {
		      console.log("couldn't handle color: " + name);
		      c = palette.black;
		      palette[name] = c;
	       }
        }
    }
    return c;
}

function makeColourSteps(steps, stops, colours) {
    var dcolours = [];
    for (var ci = 0; ci < colours.length; ++ci) {
        dcolours.push(dasColourForName(colours[ci]));
    }

    var grad = [];
  STEP_LOOP:
    for (var si = 0; si < steps; ++si) {
        var rs = (1.0 * si) / (steps-1);
        var score = stops[0] + (stops[stops.length -1] - stops[0]) * rs;
        for (var i = 0; i < stops.length - 1; ++i) {
            if (score >= stops[i] && score <= stops[i+1]) {
                var frac = (score - stops[i]) / (stops[i+1] - stops[i]);
                var ca = dcolours[i];
                var cb = dcolours[i+1];

                var fill = new DColour(
                    ((ca.red * (1.0 - frac)) + (cb.red * frac))|0,
                    ((ca.green * (1.0 - frac)) + (cb.green * frac))|0,
                    ((ca.blue * (1.0 - frac)) + (cb.blue * frac))|0
                ).toSvgString();
                grad.push(fill);

                continue STEP_LOOP;
            }
        }
        throw 'Bad step';
    }

    return grad;
}

function makeGradient(steps, color1, color2, color3) {
    if (color3) {
        return makeColourSteps(steps, [0, 0.5, 1], [color1, color2, color3]);
    } else {
        return makeColourSteps(steps, [0, 1], [color1, color2]);
    }
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        makeColourSteps: makeColourSteps,
        makeGradient: makeGradient,
        dasColourForName: dasColourForName
    };
}

},{}],5:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// das.js: queries and low-level data model.
//

"use strict";

if (typeof(require) !== 'undefined') {
    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;
    var pusho = utils.pusho;

    var color = require('./color');
    var makeColourSteps = color.makeColourSteps;
}

var dasLibErrorHandler = function(errMsg) {
    alert(errMsg);
}
var dasLibRequestQueue = new Array();



function DASSegment(name, start, end, description) {
    this.name = name;
    this.start = start;
    this.end = end;
    this.description = description;
}
DASSegment.prototype.toString = function() {
    return this.name + ':' + this.start + '..' + this.end;
};
DASSegment.prototype.isBounded = function() {
    return this.start && this.end;
}
DASSegment.prototype.toDASQuery = function() {
    var q = 'segment=' + this.name;
    if (this.start && this.end) {
        q += (':' + this.start + ',' + this.end);
    }
    return q;
}


function DASSource(a1, a2) {
    var options;
    if (typeof a1 == 'string') {
        this.uri = a1;
        options = a2 || {};
    } else {
        options = a1 || {};
    }
    for (var k in options) {
        if (typeof(options[k]) != 'function') {
            this[k] = options[k];
        }
    }


    if (!this.coords) {
        this.coords = [];
    }
    if (!this.props) {
        this.props = {};
    }

    this.dasBaseURI = this.uri;
    if (this.dasBaseURI && this.dasBaseURI.substr(this.uri.length - 1) != '/') {
        this.dasBaseURI = this.dasBaseURI + '/';
    }
}

function DASCoords() {
}

function coordsMatch(c1, c2) {
    return c1.taxon == c2.taxon && c1.auth == c2.auth && c1.version == c2.version;
}

//
// DAS 1.6 entry_points command
//

DASSource.prototype.entryPoints = function(callback) {
    var dasURI = this.dasBaseURI + 'entry_points';
    this.doCrossDomainRequest(dasURI, function(responseXML) {
            if (!responseXML) {
                return callback([]);
            }

                var entryPoints = new Array();
                
                var segs = responseXML.getElementsByTagName('SEGMENT');
                for (var i = 0; i < segs.length; ++i) {
                    var seg = segs[i];
                    var segId = seg.getAttribute('id');
                    
                    var segSize = seg.getAttribute('size');
                    var segMin, segMax;
                    if (segSize) {
                        segMin = 1; segMax = segSize|0;
                    } else {
                        segMin = seg.getAttribute('start');
                        if (segMin) {
                            segMin |= 0;
                        }
                        segMax = seg.getAttribute('stop');
                        if (segMax) {
                            segMax |= 0;
                        }
                    }
                    var segDesc = null;
                    if (seg.firstChild) {
                        segDesc = seg.firstChild.nodeValue;
                    }
                    entryPoints.push(new DASSegment(segId, segMin, segMax, segDesc));
                }          
               callback(entryPoints);
    });         
}

//
// DAS 1.6 sequence command
// Do we need an option to fall back to the dna command?
//

function DASSequence(name, start, end, alpha, seq) {
    this.name = name;
    this.start = start;
    this.end = end;
    this.alphabet = alpha;
    this.seq = seq;
}

DASSource.prototype.sequence = function(segment, callback) {
    var dasURI = this.dasBaseURI + 'sequence?' + segment.toDASQuery();
    this.doCrossDomainRequest(dasURI, function(responseXML) {
        if (!responseXML) {
            callback([]);
            return;
        } else {
                var seqs = new Array();
                
                var segs = responseXML.getElementsByTagName('SEQUENCE');
                for (var i = 0; i < segs.length; ++i) {
                    var seg = segs[i];
                    var segId = seg.getAttribute('id');
                    var segMin = seg.getAttribute('start');
                    var segMax = seg.getAttribute('stop');
                    var segAlpha = 'DNA';
                    var segSeq = null;
                    if (seg.firstChild) {
                        var rawSeq = seg.firstChild.nodeValue;
                        segSeq = '';
                        var idx = 0;
                        while (true) {
                            var space = rawSeq.indexOf('\n', idx);
                            if (space >= 0) {
                                segSeq += rawSeq.substring(idx, space).toUpperCase();
                                idx = space + 1;
                            } else {
                                segSeq += rawSeq.substring(idx).toUpperCase();
                                break;
                            }
                        }
                    }
                    seqs.push(new DASSequence(segId, segMin, segMax, segAlpha, segSeq));
                }
                
                callback(seqs);
        }
    });
}

//
// DAS 1.6 features command
//

function DASFeature() {
}

function DASGroup(id) {
    if (id)
        this.id = id;
}

function DASLink(desc, uri) {
    this.desc = desc;
    this.uri = uri;
}

DASSource.prototype.features = function(segment, options, callback) {
    options = options || {};
    var thisB = this;

    var dasURI;
    if (this.features_uri) {
        dasURI = this.features_uri;
    } else {
        var filters = [];

        if (segment) {
            filters.push(segment.toDASQuery());
        } else if (options.group) {
            var g = options.group;
            if (typeof g == 'string') {
                filters.push('group_id=' + g);
            } else {
                for (var gi = 0; gi < g.length; ++gi) {
                    filters.push('group_id=' + g[gi]);
                }
            }
        }

        if (options.adjacent) {
            var adj = options.adjacent;
            if (typeof adj == 'string') {
                adj = [adj];
            }
            for (var ai = 0; ai < adj.length; ++ai) {
                filters.push('adjacent=' + adj[ai]);
            }
        }

        if (options.type) {
            if (typeof options.type == 'string') {
                filters.push('type=' + options.type);
            } else {
                for (var ti = 0; ti < options.type.length; ++ti) {
                    filters.push('type=' + options.type[ti]);
                }
            }
        }
        
        if (options.maxbins) {
            filters.push('maxbins=' + options.maxbins);
        }
        
        if (filters.length > 0) {
            dasURI = this.dasBaseURI + 'features?' + filters.join(';');
        } else {
            callback([], 'No filters specified');
        }
    } 
   

    this.doCrossDomainRequest(dasURI, function(responseXML, req) {
        if (!responseXML) {
            var msg;
            if (req.status == 0) {
                msg = 'server may not support CORS';
            } else {
                msg = 'status=' + req.status;
            }
            callback([], 'Failed request: ' + msg);
            return;
        }
/*      if (req) {
            var caps = req.getResponseHeader('X-DAS-Capabilties');
            if (caps) {
                alert(caps);
            }
        } */

        var features = new Array();
        var segmentMap = {};

        var segs = responseXML.getElementsByTagName('SEGMENT');
        for (var si = 0; si < segs.length; ++si) {
            var segmentXML = segs[si];
            var segmentID = segmentXML.getAttribute('id');
            segmentMap[segmentID] = {
                min: segmentXML.getAttribute('start'),
                max: segmentXML.getAttribute('stop')
            };
            
            var featureXMLs = segmentXML.getElementsByTagName('FEATURE');
            for (var i = 0; i < featureXMLs.length; ++i) {
                var feature = featureXMLs[i];
                var dasFeature = new DASFeature();
                
                dasFeature.segment = segmentID;
                dasFeature.id = feature.getAttribute('id');
                dasFeature.label = feature.getAttribute('label');


/*
                var childNodes = feature.childNodes;
                for (var c = 0; c < childNodes.length; ++c) {
                    var cn = childNodes[c];
                    if (cn.nodeType == Node.ELEMENT_NODE) {
                        var key = cn.tagName;
                        //var val = null;
                        //if (cn.firstChild) {
                        //   val = cn.firstChild.nodeValue;
                        //}
                        dasFeature[key] = 'x';
                    }
                } */


                var spos = elementValue(feature, "START");
                var epos = elementValue(feature, "END");
                if ((spos|0) > (epos|0)) {
                    dasFeature.min = epos|0;
                    dasFeature.max = spos|0;
                } else {
                    dasFeature.min = spos|0;
                    dasFeature.max = epos|0;
                }
                {
                    var tec = feature.getElementsByTagName('TYPE');
                    if (tec.length > 0) {
                        var te = tec[0];
                        if (te.firstChild) {
                            dasFeature.type = te.firstChild.nodeValue;
                        }
                        dasFeature.typeId = te.getAttribute('id');
                        dasFeature.typeCv = te.getAttribute('cvId');
                    }
                }
                dasFeature.type = elementValue(feature, "TYPE");
                if (!dasFeature.type && dasFeature.typeId) {
                    dasFeature.type = dasFeature.typeId; // FIXME?
                }
                
                dasFeature.method = elementValue(feature, "METHOD");
                {
                    var ori = elementValue(feature, "ORIENTATION");
                    if (!ori) {
                        ori = '0';
                    }
                    dasFeature.orientation = ori;
                }
                dasFeature.score = elementValue(feature, "SCORE");
                dasFeature.links = dasLinksOf(feature);
                dasFeature.notes = dasNotesOf(feature);
                
                var groups = feature.getElementsByTagName("GROUP");
                for (var gi  = 0; gi < groups.length; ++gi) {
                    var groupXML = groups[gi];
                    var dasGroup = new DASGroup();
                    dasGroup.type = groupXML.getAttribute('type');
                    dasGroup.id = groupXML.getAttribute('id');
                    dasGroup.links = dasLinksOf(groupXML);
                    dasGroup.notes = dasNotesOf(groupXML);
                    if (!dasFeature.groups) {
                        dasFeature.groups = new Array(dasGroup);
                    } else {
                        dasFeature.groups.push(dasGroup);
                    }
                }

                // Magic notes.  Check with TAD before changing this.
                if (dasFeature.notes) {
                    for (var ni = 0; ni < dasFeature.notes.length; ++ni) {
                        var n = dasFeature.notes[ni];
                        if (n.indexOf('Genename=') == 0) {
                            var gg = new DASGroup();
                            gg.type='gene';
                            gg.id = n.substring(9);
                            if (!dasFeature.groups) {
                                dasFeature.groups = new Array(gg);
                            } else {
                                dasFeature.groups.push(gg);
                            }
                        }
                    }
                }
                
                {
                    var pec = feature.getElementsByTagName('PART');
                    if (pec.length > 0) {
                        var parts = [];
                        for (var pi = 0; pi < pec.length; ++pi) {
                            parts.push(pec[pi].getAttribute('id'));
                        }
                        dasFeature.parts = parts;
                    }
                }
                {
                    var pec = feature.getElementsByTagName('PARENT');
                    if (pec.length > 0) {
                        var parents = [];
                        for (var pi = 0; pi < pec.length; ++pi) {
                            parents.push(pec[pi].getAttribute('id'));
                        }
                        dasFeature.parents = parents;
                    }
                }
                
                features.push(dasFeature);
            }
        }
                
        callback(features, undefined, segmentMap);
    },
    function (err) {
        callback([], err);
    });
}

function DASAlignment(type) {
    this.type = type;
    this.objects = {};
    this.blocks = [];
}

DASSource.prototype.alignments = function(segment, options, callback) {
    var dasURI = this.dasBaseURI + 'alignment?query=' + segment;
    this.doCrossDomainRequest(dasURI, function(responseXML) {
        if (!responseXML) {
            callback([], 'Failed request ' + dasURI);
            return;
        }

        var alignments = [];
        var aliXMLs = responseXML.getElementsByTagName('alignment');
        for (var ai = 0; ai < aliXMLs.length; ++ai) {
            var aliXML = aliXMLs[ai];
            var ali = new DASAlignment(aliXML.getAttribute('alignType'));
            var objXMLs = aliXML.getElementsByTagName('alignObject');
            for (var oi = 0; oi < objXMLs.length; ++oi) {
                var objXML = objXMLs[oi];
                var obj = {
                    id:          objXML.getAttribute('intObjectId'),
                    accession:   objXML.getAttribute('dbAccessionId'),
                    version:     objXML.getAttribute('objectVersion'),
                    dbSource:    objXML.getAttribute('dbSource'),
                    dbVersion:   objXML.getAttribute('dbVersion')
                };
                ali.objects[obj.id] = obj;
            }
            
            var blockXMLs = aliXML.getElementsByTagName('block');
            for (var bi = 0; bi < blockXMLs.length; ++bi) {
                var blockXML = blockXMLs[bi];
                var block = {
                    order:      blockXML.getAttribute('blockOrder'),
                    segments:   []
                };
                var segXMLs = blockXML.getElementsByTagName('segment');
                for (var si = 0; si < segXMLs.length; ++si) {
                    var segXML = segXMLs[si];
                    var seg = {
                        object:      segXML.getAttribute('intObjectId'),
                        min:         segXML.getAttribute('start'),
                        max:         segXML.getAttribute('end'),
                        strand:      segXML.getAttribute('strand'),
                        cigar:       elementValue(segXML, 'cigar')
                    };
                    block.segments.push(seg);
                }
                ali.blocks.push(block);
            }       
                    
            alignments.push(ali);
        }
        callback(alignments);
    });
}


function DASStylesheet() {
/*
    this.highZoomStyles = new Object();
    this.mediumZoomStyles = new Object();
    this.lowZoomStyles = new Object();
*/

    this.styles = [];
}

DASStylesheet.prototype.pushStyle = function(filters, zoom, style) {
    /*

    if (!zoom) {
        this.highZoomStyles[type] = style;
        this.mediumZoomStyles[type] = style;
        this.lowZoomStyles[type] = style;
    } else if (zoom == 'high') {
        this.highZoomStyles[type] = style;
    } else if (zoom == 'medium') {
        this.mediumZoomStyles[type] = style;
    } else if (zoom == 'low') {
        this.lowZoomStyles[type] = style;
    }

    */

    if (!filters) {
        filters = {type: 'default'};
    }
    var styleHolder = shallowCopy(filters);
    if (zoom) {
        styleHolder.zoom = zoom;
    }
    styleHolder.style = style;
    this.styles.push(styleHolder);
}

function DASStyle() {
}

function parseGradient(grad) {
    var steps = grad.getAttribute('steps');
    if (steps) {
        steps = steps|0;
    } else {
        steps = 50;
    }


    var stops = [];
    var colors = [];
    var se = grad.getElementsByTagName('STOP');
    for (var si = 0; si < se.length; ++si) {
        var stop = se[si];
        stops.push(1.0 * stop.getAttribute('score'));
        colors.push(stop.firstChild.nodeValue);
    }

    return makeColourSteps(steps, stops, colors);
}

DASSource.prototype.stylesheet = function(successCB, failureCB) {
    var dasURI, creds = this.credentials;
    if (this.stylesheet_uri) {
        dasURI = this.stylesheet_uri;
        creds = false;
    } else {
        dasURI = this.dasBaseURI + 'stylesheet';
    }

    doCrossDomainRequest(dasURI, function(responseXML) {
        if (!responseXML) {
            if (failureCB) {
                failureCB();
            } 
            return;
        }
        var stylesheet = new DASStylesheet();
        var typeXMLs = responseXML.getElementsByTagName('TYPE');
        for (var i = 0; i < typeXMLs.length; ++i) {
            var typeStyle = typeXMLs[i];
            
            var filter = {};
            filter.type = typeStyle.getAttribute('id'); // Am I right in thinking that this makes DASSTYLE XML invalid?  Ugh.
            filter.label = typeStyle.getAttribute('label');
            filter.method = typeStyle.getAttribute('method');
            var glyphXMLs = typeStyle.getElementsByTagName('GLYPH');
            for (var gi = 0; gi < glyphXMLs.length; ++gi) {
                var glyphXML = glyphXMLs[gi];
                var zoom = glyphXML.getAttribute('zoom');
                var glyph = childElementOf(glyphXML);
                var style = new DASStyle();
                style.glyph = glyph.localName;
                var child = glyph.firstChild;
        
                while (child) {
                    if (child.nodeType == Node.ELEMENT_NODE) {
                        // alert(child.localName);
                        if (child.localName == 'BGGRAD') {
                            style[child.localName] = parseGradient(child);
                        } else {      
                            style[child.localName] = child.firstChild.nodeValue;
                        }
                    }
                    child = child.nextSibling;
                }
                stylesheet.pushStyle(filter, zoom, style);
            }
        }
        successCB(stylesheet);
    }, creds);
}

//
// sources command
// 

function DASRegistry(uri, opts)
{
    opts = opts || {};
    this.uri = uri;
    this.opts = opts;   
}

DASRegistry.prototype.sources = function(callback, failure, opts)
{
    if (!opts) {
        opts = {};
    }

    var filters = [];
    if (opts.taxon) {
        filters.push('organism=' + opts.taxon);
    }
    if (opts.auth) {
        filters.push('authority=' + opts.auth);
    }
    if (opts.version) {
        filters.push('version=' + opts.version);
    }
    var quri = this.uri;
    if (filters.length > 0) {
        quri = quri + '?' + filters.join('&');   // '&' as a separator to hack around dasregistry.org bug.
    }

    doCrossDomainRequest(quri, function(responseXML) {
        if (!responseXML && failure) {
            failure();
            return;
        }

        var sources = [];       
        var sourceXMLs = responseXML.getElementsByTagName('SOURCE');
        for (var si = 0; si < sourceXMLs.length; ++si) {
            var sourceXML = sourceXMLs[si];
            var versionXMLs = sourceXML.getElementsByTagName('VERSION');
            if (versionXMLs.length < 1) {
                continue;
            }
            var versionXML = versionXMLs[0];

            var coordXMLs = versionXML.getElementsByTagName('COORDINATES');
            var coords = [];
            for (var ci = 0; ci < coordXMLs.length; ++ci) {
                var coordXML = coordXMLs[ci];
                var coord = new DASCoords();
                coord.auth = coordXML.getAttribute('authority');
                coord.taxon = coordXML.getAttribute('taxid');
                coord.version = coordXML.getAttribute('version');
                coords.push(coord);
            }
            
            var caps = [];
            var capXMLs = versionXML.getElementsByTagName('CAPABILITY');
            var uri;
            for (var ci = 0; ci < capXMLs.length; ++ci) {
                var capXML = capXMLs[ci];
                
                caps.push(capXML.getAttribute('type'));

                if (capXML.getAttribute('type') == 'das1:features') {
                    var fep = capXML.getAttribute('query_uri');
                    uri = fep.substring(0, fep.length - ('features'.length));
                }
            }

            var props = {};
            var propXMLs = versionXML.getElementsByTagName('PROP');
            for (var pi = 0; pi < propXMLs.length; ++pi) {
                pusho(props, propXMLs[pi].getAttribute('name'), propXMLs[pi].getAttribute('value'));
            }
            
            if (uri) {
                var source = new DASSource(uri, {
                    source_uri: sourceXML.getAttribute('uri'),
                    name:  sourceXML.getAttribute('title'),
                    desc:  sourceXML.getAttribute('description'),
                    coords: coords,
                    props: props,
                    capabilities: caps
                });
                sources.push(source);
            }
        }
        
        callback(sources);
    });
}


//
// Utility functions
//

function elementValue(element, tag)
{
    var children = element.getElementsByTagName(tag);
    if (children.length > 0 && children[0].firstChild) {
        var c = children[0];
        if (c.childNodes.length == 1) {
            return c.firstChild.nodeValue;
        } else {
            var s = '';
            for (var ni = 0; ni < c.childNodes.length; ++ni) {
                s += c.childNodes[ni].nodeValue;
            }
            return s;
        }

    } else {
        return null;
    }
}

function childElementOf(element)
{
    if (element.hasChildNodes()) {
        var child = element.firstChild;
        do {
            if (child.nodeType == Node.ELEMENT_NODE) {
                return child;
            } 
            child = child.nextSibling;
        } while (child != null);
    }
    return null;
}


function dasLinksOf(element)
{
    var links = new Array();
    var maybeLinkChilden = element.getElementsByTagName('LINK');
    for (var ci = 0; ci < maybeLinkChilden.length; ++ci) {
        var linkXML = maybeLinkChilden[ci];
        if (linkXML.parentNode == element) {
            links.push(new DASLink(linkXML.firstChild ? linkXML.firstChild.nodeValue : 'Unknown', linkXML.getAttribute('href')));
        }
    }
    
    return links;
}

function dasNotesOf(element)
{
    var notes = [];
    var maybeNotes = element.getElementsByTagName('NOTE');
    for (var ni = 0; ni < maybeNotes.length; ++ni) {
        if (maybeNotes[ni].firstChild) {
            notes.push(maybeNotes[ni].firstChild.nodeValue);
        }
    }
    return notes;
}

function doCrossDomainRequest(url, handler, credentials, custAuth) {
    // TODO: explicit error handlers?

    if (window.XDomainRequest) {
        var req = new XDomainRequest();
        req.onload = function() {
            var dom = new ActiveXObject("Microsoft.XMLDOM");
            dom.async = false;
            dom.loadXML(req.responseText);
            handler(dom);
        }
        req.open("get", url);
        req.send('');
    } else {
        var reqStart = Date.now();
        var req = new XMLHttpRequest();

        req.onreadystatechange = function() {
            if (req.readyState == 4) {
              if (req.status >= 200 || req.status == 0) {
                  handler(req.responseXML, req);
              }
            }
        };
        req.open("get", url, true);
        if (credentials) {
            req.withCredentials = true;
        }
        if (custAuth) {
            req.setRequestHeader('X-DAS-Authorisation', custAuth);
        }
        req.overrideMimeType('text/xml');
        req.setRequestHeader('Accept', 'application/xml,*/*');
        req.send('');
    }
}

DASSource.prototype.doCrossDomainRequest = function(url, handler, errHandler) {
    var custAuth;
    if (this.xUser) {
        custAuth = 'Basic ' + btoa(this.xUser + ':' + this.xPass);
    }

    try {
        return doCrossDomainRequest(url, handler, this.credentials, custAuth);
    } catch (err) {
        if (errHandler) {
            errHandler(err);
        } else {
            throw err;
        }
    }
}

function isDasBooleanTrue(s) {
    s = ('' + s).toLowerCase();
    return s==='yes' || s==='true';
}

function isDasBooleanNotFalse(s) {
    if (!s)
        return false;

    s = ('' + s).toLowerCase();
    return s!=='no' || s!=='false';
}

function copyStylesheet(ss) {
    var nss = shallowCopy(ss);
    nss.styles = [];
    for (var si = 0; si < ss.styles.length; ++si) {
        var sh = nss.styles[si] = shallowCopy(ss.styles[si]);
        sh._methodRE = sh._labelRE = sh._typeRE = undefined;
        sh.style = shallowCopy(sh.style);
        sh.style.id = undefined;
        sh.style._gradient = undefined;
    }
    return nss;
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        DASGroup: DASGroup,
        DASFeature: DASFeature,
        DASStylesheet: DASStylesheet,
        DASStyle: DASStyle,
        DASSource: DASSource,
        DASSegment: DASSegment,
        DASRegistry: DASRegistry,
        DASSequence: DASSequence,
        DASLink: DASLink,

        isDasBooleanTrue: isDasBooleanTrue,
        isDasBooleanNotFalse: isDasBooleanNotFalse,
        copyStylesheet: copyStylesheet
    };
}
},{"./color":4,"./utils":10}],6:[function(require,module,exports){
(function (global){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// fetchworker.js
//

"use strict";

var bin = require('./bin');
var bam = require('./bam');
var bigwig = require('./bigwig');

var connections = {};

var idSeed = 0;

global.newID = function() {
    return 'cn' + (++idSeed);
}

postMessage({tag: 'init'});

self.onmessage = function(event) {
    var d = event.data;
    var command = event.data.command;
    var tag = event.data.tag;

    if (command === 'connectBAM') {
        var id = newID();

        var bamF, baiF;
        if (d.blob) {
            bamF = new bin.BlobFetchable(d.blob);
            baiF = new bin.BlobFetchable(d.indexBlob);
        } else {
            bamF = new bin.URLFetchable(d.uri, {credentials: d.credentials});
            baiF = new bin.URLFetchable(d.indexUri, {credentials: d.credentials});
        }

        bam.makeBam(bamF, baiF, function(bamObj, err) {
            if (bamObj) {
                connections[id] = new BAMWorkerFetcher(bamObj);
                postMessage({tag: tag, result: id});
            } else {
                postMessage({tag: tag, error: err || "Couldn't fetch BAM"});
            }
        });
    } else if (command === 'connectBBI') {
        var id = newID();
        var bbi;
        if (d.blob) {
            bbi = new bin.BlobFetchable(d.blob);
        } else {
            bbi = new bin.URLFetchable(d.uri, {credentials: d.credentials});
        }

        bigwig.makeBwg(bbi, function(bwg, err) {
            if (bwg) {
                connections[id] = new BBIWorkerFetcher(bwg);
                postMessage({tag: tag, result: id});
            } else {
                postMessage({tag: tag, error: err || "Couldn't fetch BBI"});
            }
        }, d.uri);
    } else if (command === 'fetch') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.fetch(d.tag, d.chr, d.min, d.max, d.zoom, d.opts);
    } else if (command === 'leap') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.leap(d.tag, d.chr, d.pos, d.dir);
    } else if (command === 'quantLeap') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.quantLeap(d.tag, d.chr, d.pos, d.dir, d.threshold, d.under);
    } else if (command === 'meta') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.meta(d.tag);
    } else if (command === 'search') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.search(d.tag, d.query, d.index);
    } else if (command === 'date') {
        return postMessage({tag: tag, result: Date.now()|0});
    } else {
        postMessage({tag: tag, error: 'Bad command ' + command});
    }
}

function BAMWorkerFetcher(bam) {
    this.bam = bam;
}

BAMWorkerFetcher.prototype.fetch = function(tag, chr, min, max, zoom, opts) {
    opts = opts || {};
    this.bam.fetch(chr, min, max, function(records, err) {
        if (records) {
            postMessage({tag: tag, result: records, time: Date.now()|0});
        } else {
            postMessage({tag: tag, error: err});
        }
    }, opts);
}

function BBIWorkerFetcher(bbi) {
    this.bbi = bbi;
}

BBIWorkerFetcher.prototype.fetch = function(tag, chr, min, max, zoom) {
    if (typeof(zoom) !== 'number')
        zoom = -1;

    var data;
    if (zoom < 0) {
        data = this.bbi.getUnzoomedView();
    } else {
        data = this.bbi.getZoomedView(zoom);
    }

    data.readWigData(chr, min, max, function(features) {
        postMessage({tag: tag, result: features});
    });
}

BBIWorkerFetcher.prototype.meta = function(tag) {
    var scales = [1];
    for (var z = 0; z < this.bbi.zoomLevels.length; ++z) {
        scales.push(this.bbi.zoomLevels[z].reduction);
    }

    var thisB = this;
    var meta = {type: this.bbi.type,
                zoomLevels: scales,
                fieldCount: this.bbi.fieldCount,
                definedFieldCount: this.bbi.definedFieldCount,
                schema: this.bbi.schema};
    if (this.bbi.type === 'bigbed') {
        this.bbi.getExtraIndices(function(ei) {
            if (ei) {
                thisB.extraIndices = ei;
                meta.extraIndices = ei.map(function(i) {return i.field});
            }
            postMessage({tag: tag, result: meta});
        });
    } else {
        postMessage({tag: tag, result: meta});
    }
}

BBIWorkerFetcher.prototype.leap = function(tag, chr, pos, dir) {
    this.bbi.getUnzoomedView().getFirstAdjacent(chr, pos, dir, function(result, err) {
        postMessage({tag: tag, result: result, error: err});
    });
}

BBIWorkerFetcher.prototype.quantLeap = function(tag, chr, pos, dir, threshold, under) {
    this.bbi.thresholdSearch(chr, pos, dir, threshold, function(result, err) {
        postMessage({tag: tag, result: result, error: err});
    });
}

BBIWorkerFetcher.prototype.search = function(tag, query, index) {
    var is = this.extraIndices[0];
    is.lookup(query, function(result, err) {
        postMessage({tag: tag, result: result, error: err});
    });
}

}).call(this,typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./bam":1,"./bigwig":2,"./bin":3}],7:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// lh3utils.js: common support for lh3's file formats
//

if (typeof(require) !== 'undefined') {
    var jszlib = require('jszlib');
    var jszlib_inflate_buffer = jszlib.inflateBuffer;
    var arrayCopy = jszlib.arrayCopy;
}

function Vob(b, o) {
    this.block = b;
    this.offset = o;
}

Vob.prototype.toString = function() {
    return '' + this.block + ':' + this.offset;
}

function readVob(ba, offset) {
    var block = ((ba[offset+6] & 0xff) * 0x100000000) + ((ba[offset+5] & 0xff) * 0x1000000) + ((ba[offset+4] & 0xff) * 0x10000) + ((ba[offset+3] & 0xff) * 0x100) + ((ba[offset+2] & 0xff));
    var bint = (ba[offset+1] << 8) | (ba[offset]);
    if (block == 0 && bint == 0) {
        return null;  // Should only happen in the linear index?
    } else {
        return new Vob(block, bint);
    }
}

function unbgzf(data, lim) {
    lim = Math.min(lim || 1, data.byteLength - 50);
    var oBlockList = [];
    var ptr = [0];
    var totalSize = 0;

    while (ptr[0] < lim) {
        var ba = new Uint8Array(data, ptr[0], 12); // FIXME is this enough for all credible BGZF block headers?
        var xlen = (ba[11] << 8) | (ba[10]);
        // dlog('xlen[' + (ptr[0]) +']=' + xlen);
        var unc = jszlib_inflate_buffer(data, 12 + xlen + ptr[0], Math.min(65536, data.byteLength - 12 - xlen - ptr[0]), ptr);
        ptr[0] += 8;
        totalSize += unc.byteLength;
        oBlockList.push(unc);
    }

    if (oBlockList.length == 1) {
        return oBlockList[0];
    } else {
        var out = new Uint8Array(totalSize);
        var cursor = 0;
        for (var i = 0; i < oBlockList.length; ++i) {
            var b = new Uint8Array(oBlockList[i]);
            arrayCopy(b, 0, out, cursor, b.length);
            cursor += b.length;
        }
        return out.buffer;
    }
}

function Chunk(minv, maxv) {
    this.minv = minv; this.maxv = maxv;
}


//
// Binning (transliterated from SAM1.3 spec)
//

/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
function reg2bin(beg, end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
var MAX_BIN = (((1<<18)-1)/7);
function reg2bins(beg, end) 
{
    var i = 0, k, list = [];
    --end;
    list.push(0);
    for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) list.push(k);
    for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) list.push(k);
    for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) list.push(k);
    for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) list.push(k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list.push(k);
    return list;
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        unbgzf: unbgzf,
        readVob: readVob,
        reg2bin: reg2bin,
        reg2bins: reg2bins,
        Chunk: Chunk
    };
}
},{"jszlib":11}],8:[function(require,module,exports){
/*
 * A JavaScript implementation of the Secure Hash Algorithm, SHA-1, as defined
 * in FIPS 180-1
 * Version 2.2 Copyright Paul Johnston 2000 - 2009.
 * Other contributors: Greg Holt, Andrew Kepert, Ydnar, Lostinet
 * Distributed under the BSD License
 * See http://pajhome.org.uk/crypt/md5 for details.
 */

 "use strict";

/*
 * Configurable variables. You may need to tweak these to be compatible with
 * the server-side, but the defaults work in most cases.
 */
var hexcase = 0;  /* hex output format. 0 - lowercase; 1 - uppercase        */
var b64pad  = ""; /* base-64 pad character. "=" for strict RFC compliance   */

/*
 * These are the functions you'll usually want to call
 * They take string arguments and return either hex or base-64 encoded strings
 */
function hex_sha1(s)    { return rstr2hex(rstr_sha1(str2rstr_utf8(s))); }
function b64_sha1(s)    { return rstr2b64(rstr_sha1(str2rstr_utf8(s))); }
function any_sha1(s, e) { return rstr2any(rstr_sha1(str2rstr_utf8(s)), e); }
function hex_hmac_sha1(k, d)
  { return rstr2hex(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d))); }
function b64_hmac_sha1(k, d)
  { return rstr2b64(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d))); }
function any_hmac_sha1(k, d, e)
  { return rstr2any(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d)), e); }

/*
 * Perform a simple self-test to see if the VM is working
 */
function sha1_vm_test()
{
  return hex_sha1("abc").toLowerCase() == "a9993e364706816aba3e25717850c26c9cd0d89d";
}

/*
 * Calculate the SHA1 of a raw string
 */
function rstr_sha1(s)
{
  return binb2rstr(binb_sha1(rstr2binb(s), s.length * 8));
}

/*
 * Calculate the HMAC-SHA1 of a key and some data (raw strings)
 */
function rstr_hmac_sha1(key, data)
{
  var bkey = rstr2binb(key);
  if(bkey.length > 16) bkey = binb_sha1(bkey, key.length * 8);

  var ipad = Array(16), opad = Array(16);
  for(var i = 0; i < 16; i++)
  {
    ipad[i] = bkey[i] ^ 0x36363636;
    opad[i] = bkey[i] ^ 0x5C5C5C5C;
  }

  var hash = binb_sha1(ipad.concat(rstr2binb(data)), 512 + data.length * 8);
  return binb2rstr(binb_sha1(opad.concat(hash), 512 + 160));
}

/*
 * Convert a raw string to a hex string
 */
function rstr2hex(input)
{
  // try { hexcase } catch(e) { hexcase=0; }
  var hex_tab = hexcase ? "0123456789ABCDEF" : "0123456789abcdef";
  var output = "";
  var x;
  for(var i = 0; i < input.length; i++)
  {
    x = input.charCodeAt(i);
    output += hex_tab.charAt((x >>> 4) & 0x0F)
           +  hex_tab.charAt( x        & 0x0F);
  }
  return output;
}

/*
 * Convert a raw string to a base-64 string
 */
function rstr2b64(input)
{
  // try { b64pad } catch(e) { b64pad=''; }
  var tab = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  var output = "";
  var len = input.length;
  for(var i = 0; i < len; i += 3)
  {
    var triplet = (input.charCodeAt(i) << 16)
                | (i + 1 < len ? input.charCodeAt(i+1) << 8 : 0)
                | (i + 2 < len ? input.charCodeAt(i+2)      : 0);
    for(var j = 0; j < 4; j++)
    {
      if(i * 8 + j * 6 > input.length * 8) output += b64pad;
      else output += tab.charAt((triplet >>> 6*(3-j)) & 0x3F);
    }
  }
  return output;
}

/*
 * Convert a raw string to an arbitrary string encoding
 */
function rstr2any(input, encoding)
{
  var divisor = encoding.length;
  var remainders = Array();
  var i, q, x, quotient;

  /* Convert to an array of 16-bit big-endian values, forming the dividend */
  var dividend = Array(Math.ceil(input.length / 2));
  for(i = 0; i < dividend.length; i++)
  {
    dividend[i] = (input.charCodeAt(i * 2) << 8) | input.charCodeAt(i * 2 + 1);
  }

  /*
   * Repeatedly perform a long division. The binary array forms the dividend,
   * the length of the encoding is the divisor. Once computed, the quotient
   * forms the dividend for the next step. We stop when the dividend is zero.
   * All remainders are stored for later use.
   */
  while(dividend.length > 0)
  {
    quotient = Array();
    x = 0;
    for(i = 0; i < dividend.length; i++)
    {
      x = (x << 16) + dividend[i];
      q = Math.floor(x / divisor);
      x -= q * divisor;
      if(quotient.length > 0 || q > 0)
        quotient[quotient.length] = q;
    }
    remainders[remainders.length] = x;
    dividend = quotient;
  }

  /* Convert the remainders to the output string */
  var output = "";
  for(i = remainders.length - 1; i >= 0; i--)
    output += encoding.charAt(remainders[i]);

  /* Append leading zero equivalents */
  var full_length = Math.ceil(input.length * 8 /
                                    (Math.log(encoding.length) / Math.log(2)))
  for(i = output.length; i < full_length; i++)
    output = encoding[0] + output;

  return output;
}

/*
 * Encode a string as utf-8.
 * For efficiency, this assumes the input is valid utf-16.
 */
function str2rstr_utf8(input)
{
  var output = "";
  var i = -1;
  var x, y;

  while(++i < input.length)
  {
    /* Decode utf-16 surrogate pairs */
    x = input.charCodeAt(i);
    y = i + 1 < input.length ? input.charCodeAt(i + 1) : 0;
    if(0xD800 <= x && x <= 0xDBFF && 0xDC00 <= y && y <= 0xDFFF)
    {
      x = 0x10000 + ((x & 0x03FF) << 10) + (y & 0x03FF);
      i++;
    }

    /* Encode output as utf-8 */
    if(x <= 0x7F)
      output += String.fromCharCode(x);
    else if(x <= 0x7FF)
      output += String.fromCharCode(0xC0 | ((x >>> 6 ) & 0x1F),
                                    0x80 | ( x         & 0x3F));
    else if(x <= 0xFFFF)
      output += String.fromCharCode(0xE0 | ((x >>> 12) & 0x0F),
                                    0x80 | ((x >>> 6 ) & 0x3F),
                                    0x80 | ( x         & 0x3F));
    else if(x <= 0x1FFFFF)
      output += String.fromCharCode(0xF0 | ((x >>> 18) & 0x07),
                                    0x80 | ((x >>> 12) & 0x3F),
                                    0x80 | ((x >>> 6 ) & 0x3F),
                                    0x80 | ( x         & 0x3F));
  }
  return output;
}

/*
 * Encode a string as utf-16
 */
function str2rstr_utf16le(input)
{
  var output = "";
  for(var i = 0; i < input.length; i++)
    output += String.fromCharCode( input.charCodeAt(i)        & 0xFF,
                                  (input.charCodeAt(i) >>> 8) & 0xFF);
  return output;
}

function str2rstr_utf16be(input)
{
  var output = "";
  for(var i = 0; i < input.length; i++)
    output += String.fromCharCode((input.charCodeAt(i) >>> 8) & 0xFF,
                                   input.charCodeAt(i)        & 0xFF);
  return output;
}

/*
 * Convert a raw string to an array of big-endian words
 * Characters >255 have their high-byte silently ignored.
 */
function rstr2binb(input)
{
  var output = Array(input.length >> 2);
  for(var i = 0; i < output.length; i++)
    output[i] = 0;
  for(var i = 0; i < input.length * 8; i += 8)
    output[i>>5] |= (input.charCodeAt(i / 8) & 0xFF) << (24 - i % 32);
  return output;
}

/*
 * Convert an array of big-endian words to a string
 */
function binb2rstr(input)
{
  var output = "";
  for(var i = 0; i < input.length * 32; i += 8)
    output += String.fromCharCode((input[i>>5] >>> (24 - i % 32)) & 0xFF);
  return output;
}

/*
 * Calculate the SHA-1 of an array of big-endian words, and a bit length
 */
function binb_sha1(x, len)
{
  /* append padding */
  x[len >> 5] |= 0x80 << (24 - len % 32);
  x[((len + 64 >> 9) << 4) + 15] = len;

  var w = Array(80);
  var a =  1732584193;
  var b = -271733879;
  var c = -1732584194;
  var d =  271733878;
  var e = -1009589776;

  for(var i = 0; i < x.length; i += 16)
  {
    var olda = a;
    var oldb = b;
    var oldc = c;
    var oldd = d;
    var olde = e;

    for(var j = 0; j < 80; j++)
    {
      if(j < 16) w[j] = x[i + j];
      else w[j] = bit_rol(w[j-3] ^ w[j-8] ^ w[j-14] ^ w[j-16], 1);
      var t = safe_add(safe_add(bit_rol(a, 5), sha1_ft(j, b, c, d)),
                       safe_add(safe_add(e, w[j]), sha1_kt(j)));
      e = d;
      d = c;
      c = bit_rol(b, 30);
      b = a;
      a = t;
    }

    a = safe_add(a, olda);
    b = safe_add(b, oldb);
    c = safe_add(c, oldc);
    d = safe_add(d, oldd);
    e = safe_add(e, olde);
  }
  return Array(a, b, c, d, e);

}

/*
 * Perform the appropriate triplet combination function for the current
 * iteration
 */
function sha1_ft(t, b, c, d)
{
  if(t < 20) return (b & c) | ((~b) & d);
  if(t < 40) return b ^ c ^ d;
  if(t < 60) return (b & c) | (b & d) | (c & d);
  return b ^ c ^ d;
}

/*
 * Determine the appropriate additive constant for the current iteration
 */
function sha1_kt(t)
{
  return (t < 20) ?  1518500249 : (t < 40) ?  1859775393 :
         (t < 60) ? -1894007588 : -899497514;
}

/*
 * Add integers, wrapping at 2^32. This uses 16-bit operations internally
 * to work around bugs in some JS interpreters.
 */
function safe_add(x, y)
{
  var lsw = (x & 0xFFFF) + (y & 0xFFFF);
  var msw = (x >> 16) + (y >> 16) + (lsw >> 16);
  return (msw << 16) | (lsw & 0xFFFF);
}

/*
 * Bitwise rotate a 32-bit number to the left.
 */
function bit_rol(num, cnt)
{
  return (num << cnt) | (num >>> (32 - cnt));
}

if (typeof(module) !== 'undefined') {
  module.exports = {
    b64_sha1: b64_sha1,
    hex_sha1: hex_sha1
  }
}

},{}],9:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// spans.js: JavaScript Intset/Location port.
//

"use strict";


function Range(min, max)
{
    if (typeof(min) != 'number' || typeof(max) != 'number')
        throw 'Bad range ' + min + ',' + max;
    this._min = min;
    this._max = max;
}

Range.prototype.min = function() {
    return this._min;
}

Range.prototype.max = function() {
    return this._max;
}

Range.prototype.contains = function(pos) {
    return pos >= this._min && pos <= this._max;
}

Range.prototype.isContiguous = function() {
    return true;
}

Range.prototype.ranges = function() {
    return [this];
}

Range.prototype._pushRanges = function(ranges) {
    ranges.push(this);
}

Range.prototype.toString = function() {
    return '[' + this._min + '-' + this._max + ']';
}

function _Compound(ranges) {
    this._ranges = ranges;
    // assert sorted?
}

_Compound.prototype.min = function() {
    return this._ranges[0].min();
}

_Compound.prototype.max = function() {
    return this._ranges[this._ranges.length - 1].max();
}

_Compound.prototype.contains = function(pos) {
    // FIXME implement bsearch if we use this much.
    for (var s = 0; s < this._ranges.length; ++s) {
        if (this._ranges[s].contains(pos)) {
            return true;
        }
    }
    return false;
}

_Compound.prototype.isContiguous = function() {
    return this._ranges.length > 1;
}

_Compound.prototype.ranges = function() {
    return this._ranges;
}

_Compound.prototype._pushRanges = function(ranges) {
    for (var ri = 0; ri < this._ranges.length; ++ri)
        ranges.push(this._ranges[ri]);
}

_Compound.prototype.toString = function() {
    var s = '';
    for (var r = 0; r < this._ranges.length; ++r) {
        if (r>0) {
            s = s + ',';
        }
        s = s + this._ranges[r].toString();
    }
    return s;
}

function union(s0, s1) {
    if (! (s0 instanceof Array)) {
        s0 = [s0];
        if (s1)
            s0.push(s1);
    }

    if (s0.length == 0)
        return null;
    else if (s0.length == 1)
        return s0[0];

    var ranges = [];
    for (var si = 0; si < s0.length; ++si)
        s0[si]._pushRanges(ranges);
    ranges = ranges.sort(_rangeOrder);

    var oranges = [];
    var current = ranges[0];
    current = new Range(current._min, current._max);  // Copy now so we don't have to later.

    for (var i = 1; i < ranges.length; ++i) {
        var nxt = ranges[i];
        if (nxt._min > (current._max + 1)) {
            oranges.push(current);
            current = new Range(nxt._min, nxt._max);
        } else {
            if (nxt._max > current._max) {
                current._max = nxt._max;
            }
        }
    }
    oranges.push(current);

    if (oranges.length == 1) {
        return oranges[0];
    } else {
        return new _Compound(oranges);
    }
}

function intersection(s0, s1) {
    var r0 = s0.ranges();
    var r1 = s1.ranges();
    var l0 = r0.length, l1 = r1.length;
    var i0 = 0, i1 = 0;
    var or = [];

    while (i0 < l0 && i1 < l1) {
        var s0 = r0[i0], s1 = r1[i1];
        var lapMin = Math.max(s0.min(), s1.min());
        var lapMax = Math.min(s0.max(), s1.max());
        if (lapMax >= lapMin) {
            or.push(new Range(lapMin, lapMax));
        }
        if (s0.max() > s1.max()) {
            ++i1;
        } else {
            ++i0;
        }
    }
    
    if (or.length == 0) {
        return null; // FIXME
    } else if (or.length == 1) {
        return or[0];
    } else {
        return new _Compound(or);
    }
}

function coverage(s) {
    var tot = 0;
    var rl = s.ranges();
    for (var ri = 0; ri < rl.length; ++ri) {
        var r = rl[ri];
        tot += (r.max() - r.min() + 1);
    }
    return tot;
}



function rangeOrder(a, b)
{
    if (a.min() < b.min()) {
        return -1;
    } else if (a.min() > b.min()) {
        return 1;
    } else if (a.max() < b.max()) {
        return -1;
    } else if (b.max() > a.max()) {
        return 1;
    } else {
        return 0;
    }
}

function _rangeOrder(a, b)
{
    if (a._min < b._min) {
        return -1;
    } else if (a._min > b._min) {
        return 1;
    } else if (a._max < b._max) {
        return -1;
    } else if (b._max > a._max) {
        return 1;
    } else {
        return 0;
    }
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        Range: Range,
        union: union,
        intersection: intersection,
        coverage: coverage,
        rangeOver: rangeOrder,
        _rangeOrder: _rangeOrder
    }
}
},{}],10:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// utils.js: odds, sods, and ends.
//

"use strict";

if (typeof(require) !== 'undefined') {
    var sha1 = require('./sha1');
    var b64_sha1 = sha1.b64_sha1;
}

var NUM_REGEXP = new RegExp('[0-9]+');

function stringToNumbersArray(str) {
    var nums = new Array();
    var m;
    while (m = NUM_REGEXP.exec(str)) {
        nums.push(m[0]);
        str=str.substring(m.index + (m[0].length));
    }
    return nums;
}

var STRICT_NUM_REGEXP = new RegExp('^[0-9]+$');

function stringToInt(str) {
    str = str.replace(new RegExp(',', 'g'), '');
    if (!STRICT_NUM_REGEXP.test(str)) {
        return null;
    }
    return str|0;
}

function pushnew(a, v) {
    for (var i = 0; i < a.length; ++i) {
        if (a[i] == v) {
            return;
        }
    }
    a.push(v);
}

function pusho(obj, k, v) {
    if (obj[k]) {
        obj[k].push(v);
    } else {
        obj[k] = [v];
    }
}

function pushnewo(obj, k, v) {
    var a = obj[k];
    if (a) {
        for (var i = 0; i < a.length; ++i) {    // indexOf requires JS16 :-(.
            if (a[i] == v) {
                return;
            }
        }
        a.push(v);
    } else {
        obj[k] = [v];
    }
}


function pick(a, b, c, d)
{
    if (a) {
        return a;
    } else if (b) {
        return b;
    } else if (c) {
        return c;
    } else if (d) {
        return d;
    }
}

function pushnew(l, o)
{
    for (var i = 0; i < l.length; ++i) {
        if (l[i] == o) {
            return;
        }
    }
    l.push(o);
}



function arrayIndexOf(a, x) {
    if (!a) {
        return -1;
    }

    for (var i = 0; i < a.length; ++i) {
        if (a[i] === x) {
            return i;
        }
    }
    return -1;
}

function arrayRemove(a, x) {
    var i = arrayIndexOf(a, x);
    if (i >= 0) {
        a.splice(i, 1);
        return true;
    }
    return false;
}

//
// DOM utilities
//


function makeElement(tag, children, attribs, styles)
{
    var ele = document.createElement(tag);
    if (children) {
        if (! (children instanceof Array)) {
            children = [children];
        }
        for (var i = 0; i < children.length; ++i) {
            var c = children[i];
            if (c) {
                if (typeof c == 'string') {
                    c = document.createTextNode(c);
                } else if (typeof c == 'number') {
                    c = document.createTextNode('' + c);
                }
                ele.appendChild(c);
            }
        }
    }
    
    if (attribs) {
        for (var l in attribs) {
            try {
                ele[l] = attribs[l];
            } catch (e) {
                console.log('error setting ' + l);
                throw(e);
            }
        }
    }
    if (styles) {
        for (var l in styles) {
            ele.style[l] = styles[l];
        }
    }
    return ele;
}

function makeElementNS(namespace, tag, children, attribs)
{
    var ele = document.createElementNS(namespace, tag);
    if (children) {
        if (! (children instanceof Array)) {
            children = [children];
        }
        for (var i = 0; i < children.length; ++i) {
            var c = children[i];
            if (typeof c == 'string') {
                c = document.createTextNode(c);
            }
            ele.appendChild(c);
        }
    }
    
    setAttrs(ele, attribs);
    return ele;
}

var attr_name_cache = {};

function setAttr(node, key, value)
{
    var attr = attr_name_cache[key];
    if (!attr) {
        var _attr = '';
        for (var c = 0; c < key.length; ++c) {
            var cc = key.substring(c, c+1);
            var lcc = cc.toLowerCase();
            if (lcc != cc) {
                _attr = _attr + '-' + lcc;
            } else {
                _attr = _attr + cc;
            }
        }
        attr_name_cache[key] = _attr;
        attr = _attr;
    }
    node.setAttribute(attr, value);
}

function setAttrs(node, attribs)
{
    if (attribs) {
        for (var l in attribs) {
            setAttr(node, l, attribs[l]);
        }
    }
}



function removeChildren(node)
{
    if (!node || !node.childNodes) {
        return;
    }

    while (node.childNodes.length > 0) {
        node.removeChild(node.firstChild);
    }
}



//
// WARNING: not for general use!
//

function miniJSONify(o, exc) {
    if (typeof o === 'undefined') {
        return 'undefined';
    } else if (o == null) {
        return 'null';
    } else if (typeof o == 'string') {
        return "'" + o + "'";
    } else if (typeof o == 'number') {
        return "" + o;
    } else if (typeof o == 'boolean') {
        return "" + o;
    } else if (typeof o == 'object') {
        if (o instanceof Array) {
            var s = null;
            for (var i = 0; i < o.length; ++i) {
                s = (s == null ? '' : (s + ', ')) + miniJSONify(o[i], exc);
            }
            return '[' + (s?s:'') + ']';
        } else {
            exc = exc || {};
            var s = null;
            for (var k in o) {
                if (exc[k])
                    continue;
                if (k != undefined && typeof(o[k]) != 'function') {
                    s = (s == null ? '' : (s + ', ')) + k + ': ' + miniJSONify(o[k], exc);
                }
            }
            return '{' + (s?s:'') + '}';
        }
    } else {
        return (typeof o);
    }
}

function shallowCopy(o) {
    var n = {};
    for (var k in o) {
        n[k] = o[k];
    }
    return n;
}

function Observed(x) {
    this.value = x;
    this.listeners = [];
}

Observed.prototype.addListener = function(f) {
    this.listeners.push(f);
}

Observed.prototype.addListenerAndFire = function(f) {
    this.listeners.push(f);
    f(this.value);
}

Observed.prototype.removeListener = function(f) {
    arrayRemove(this.listeners, f);
}

Observed.prototype.get = function() {
    return this.value;
}

Observed.prototype.set = function(x) {
    this.value = x;
    for (var i = 0; i < this.listeners.length; ++i) {
        this.listeners[i](x);
    }
}

function Awaited() {
    this.queue = [];
}

Awaited.prototype.provide = function(x) {
    if (this.res !== undefined) {
        throw "Resource has already been provided.";
    }

    this.res = x;
    for (var i = 0; i < this.queue.length; ++i) {
        this.queue[i](x);
    }
    this.queue = null;   // avoid leaking closures.
}

Awaited.prototype.await = function(f) {
    if (this.res !== undefined) {
        f(this.res);
        return this.res;
    } else {
        this.queue.push(f);
    }
}

var __dalliance_saltSeed = 0;

function saltURL(url) {
    return url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++__dalliance_saltSeed));
}

function textXHR(url, callback, opts) {
    if (opts.salt) 
        url = saltURL(url);

    var req = new XMLHttpRequest();
    req.onreadystatechange = function() {
    	if (req.readyState == 4) {
    	    if (req.status >= 300) {
    		    callback(null, 'Error code ' + req.status);
    	    } else {
    		    callback(req.responseText);
    	    }
    	}
    };
    
    req.open('GET', url, true);
    req.responseType = 'text';

    if (opts && opts.credentials) {
        req.withCredentials = true;
    }
    req.send('');
}

function relativeURL(base, rel) {
    // FIXME quite naive -- good enough for trackhubs?

    if (rel.indexOf('http:') == 0 || rel.indexOf('https:') == 0) {
        return rel;
    }

    var li = base.lastIndexOf('/');
    if (li >= 0) {
        return base.substr(0, li + 1) + rel;
    } else {
        return rel;
    }
}

var AMINO_ACID_TRANSLATION = {
    'TTT': 'F',
    'TTC': 'F',
    'TTA': 'L',
    'TTG': 'L',
    'CTT': 'L',
    'CTC': 'L',
    'CTA': 'L',
    'CTG': 'L',
    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'ATG': 'M',
    'GTT': 'V',
    'GTC': 'V',
    'GTA': 'V',
    'GTG': 'V',
    'TCT': 'S',
    'TCC': 'S',
    'TCA': 'S',
    'TCG': 'S',
    'CCT': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',
    'ACT': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',
    'GCT': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',
    'TAT': 'Y',
    'TAC': 'Y',
    'TAA': '*',  // stop
    'TAG': '*',  // stop
    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'AAT': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'TGT': 'C',
    'TGC': 'C',
    'TGA': '*',  // stop
    'TGG': 'W',
    'CGT': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GGT': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G'
}

function resolveUrlToPage(rel) {
    return makeElement('a', null, {href: rel}).href;
}

//
// Missing APIs
// 

if (!('trim' in String.prototype)) {
    String.prototype.trim = function() {
        return this.replace(/^\s+/, '').replace(/\s+$/, '');
    };
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        textXHR: textXHR,
        relativeURL: relativeURL,
        resolveUrlToPage: resolveUrlToPage,
        shallowCopy: shallowCopy,
        pusho: pusho,
        pushnew: pushnew,
        pushnewo: pushnewo,
        arrayIndexOf: arrayIndexOf,
        pick: pick,

        makeElement: makeElement,
        makeElementNS: makeElementNS,
        removeChildren: removeChildren,

        miniJSONify: miniJSONify,

        Observed: Observed,
        Awaited: Awaited,

        AMINO_ACID_TRANSLATION: AMINO_ACID_TRANSLATION
    }
}

},{"./sha1":8}],11:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Javascript ZLib
// By Thomas Down 2010-2011
//
// Based very heavily on portions of jzlib (by ymnk@jcraft.com), who in
// turn credits Jean-loup Gailly and Mark Adler for the original zlib code.
//
// inflate.js: ZLib inflate code
//

//
// Shared constants
//

var MAX_WBITS=15; // 32K LZ77 window
var DEF_WBITS=MAX_WBITS;
var MAX_MEM_LEVEL=9;
var MANY=1440;
var BMAX = 15;

// preset dictionary flag in zlib header
var PRESET_DICT=0x20;

var Z_NO_FLUSH=0;
var Z_PARTIAL_FLUSH=1;
var Z_SYNC_FLUSH=2;
var Z_FULL_FLUSH=3;
var Z_FINISH=4;

var Z_DEFLATED=8;

var Z_OK=0;
var Z_STREAM_END=1;
var Z_NEED_DICT=2;
var Z_ERRNO=-1;
var Z_STREAM_ERROR=-2;
var Z_DATA_ERROR=-3;
var Z_MEM_ERROR=-4;
var Z_BUF_ERROR=-5;
var Z_VERSION_ERROR=-6;

var METHOD=0;   // waiting for method byte
var FLAG=1;     // waiting for flag byte
var DICT4=2;    // four dictionary check bytes to go
var DICT3=3;    // three dictionary check bytes to go
var DICT2=4;    // two dictionary check bytes to go
var DICT1=5;    // one dictionary check byte to go
var DICT0=6;    // waiting for inflateSetDictionary
var BLOCKS=7;   // decompressing blocks
var CHECK4=8;   // four check bytes to go
var CHECK3=9;   // three check bytes to go
var CHECK2=10;  // two check bytes to go
var CHECK1=11;  // one check byte to go
var DONE=12;    // finished check, done
var BAD=13;     // got an error--stay here

var inflate_mask = [0x00000000, 0x00000001, 0x00000003, 0x00000007, 0x0000000f, 0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff, 0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff, 0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff];

var IB_TYPE=0;  // get type bits (3, including end bit)
var IB_LENS=1;  // get lengths for stored
var IB_STORED=2;// processing stored block
var IB_TABLE=3; // get table lengths
var IB_BTREE=4; // get bit lengths tree for a dynamic block
var IB_DTREE=5; // get length, distance trees for a dynamic block
var IB_CODES=6; // processing fixed or dynamic block
var IB_DRY=7;   // output remaining window bytes
var IB_DONE=8;  // finished last block, done
var IB_BAD=9;   // ot a data error--stuck here

var fixed_bl = 9;
var fixed_bd = 5;

var fixed_tl = [
    96,7,256, 0,8,80, 0,8,16, 84,8,115,
    82,7,31, 0,8,112, 0,8,48, 0,9,192,
    80,7,10, 0,8,96, 0,8,32, 0,9,160,
    0,8,0, 0,8,128, 0,8,64, 0,9,224,
    80,7,6, 0,8,88, 0,8,24, 0,9,144,
    83,7,59, 0,8,120, 0,8,56, 0,9,208,
    81,7,17, 0,8,104, 0,8,40, 0,9,176,
    0,8,8, 0,8,136, 0,8,72, 0,9,240,
    80,7,4, 0,8,84, 0,8,20, 85,8,227,
    83,7,43, 0,8,116, 0,8,52, 0,9,200,
    81,7,13, 0,8,100, 0,8,36, 0,9,168,
    0,8,4, 0,8,132, 0,8,68, 0,9,232,
    80,7,8, 0,8,92, 0,8,28, 0,9,152,
    84,7,83, 0,8,124, 0,8,60, 0,9,216,
    82,7,23, 0,8,108, 0,8,44, 0,9,184,
    0,8,12, 0,8,140, 0,8,76, 0,9,248,
    80,7,3, 0,8,82, 0,8,18, 85,8,163,
    83,7,35, 0,8,114, 0,8,50, 0,9,196,
    81,7,11, 0,8,98, 0,8,34, 0,9,164,
    0,8,2, 0,8,130, 0,8,66, 0,9,228,
    80,7,7, 0,8,90, 0,8,26, 0,9,148,
    84,7,67, 0,8,122, 0,8,58, 0,9,212,
    82,7,19, 0,8,106, 0,8,42, 0,9,180,
    0,8,10, 0,8,138, 0,8,74, 0,9,244,
    80,7,5, 0,8,86, 0,8,22, 192,8,0,
    83,7,51, 0,8,118, 0,8,54, 0,9,204,
    81,7,15, 0,8,102, 0,8,38, 0,9,172,
    0,8,6, 0,8,134, 0,8,70, 0,9,236,
    80,7,9, 0,8,94, 0,8,30, 0,9,156,
    84,7,99, 0,8,126, 0,8,62, 0,9,220,
    82,7,27, 0,8,110, 0,8,46, 0,9,188,
    0,8,14, 0,8,142, 0,8,78, 0,9,252,
    96,7,256, 0,8,81, 0,8,17, 85,8,131,
    82,7,31, 0,8,113, 0,8,49, 0,9,194,
    80,7,10, 0,8,97, 0,8,33, 0,9,162,
    0,8,1, 0,8,129, 0,8,65, 0,9,226,
    80,7,6, 0,8,89, 0,8,25, 0,9,146,
    83,7,59, 0,8,121, 0,8,57, 0,9,210,
    81,7,17, 0,8,105, 0,8,41, 0,9,178,
    0,8,9, 0,8,137, 0,8,73, 0,9,242,
    80,7,4, 0,8,85, 0,8,21, 80,8,258,
    83,7,43, 0,8,117, 0,8,53, 0,9,202,
    81,7,13, 0,8,101, 0,8,37, 0,9,170,
    0,8,5, 0,8,133, 0,8,69, 0,9,234,
    80,7,8, 0,8,93, 0,8,29, 0,9,154,
    84,7,83, 0,8,125, 0,8,61, 0,9,218,
    82,7,23, 0,8,109, 0,8,45, 0,9,186,
    0,8,13, 0,8,141, 0,8,77, 0,9,250,
    80,7,3, 0,8,83, 0,8,19, 85,8,195,
    83,7,35, 0,8,115, 0,8,51, 0,9,198,
    81,7,11, 0,8,99, 0,8,35, 0,9,166,
    0,8,3, 0,8,131, 0,8,67, 0,9,230,
    80,7,7, 0,8,91, 0,8,27, 0,9,150,
    84,7,67, 0,8,123, 0,8,59, 0,9,214,
    82,7,19, 0,8,107, 0,8,43, 0,9,182,
    0,8,11, 0,8,139, 0,8,75, 0,9,246,
    80,7,5, 0,8,87, 0,8,23, 192,8,0,
    83,7,51, 0,8,119, 0,8,55, 0,9,206,
    81,7,15, 0,8,103, 0,8,39, 0,9,174,
    0,8,7, 0,8,135, 0,8,71, 0,9,238,
    80,7,9, 0,8,95, 0,8,31, 0,9,158,
    84,7,99, 0,8,127, 0,8,63, 0,9,222,
    82,7,27, 0,8,111, 0,8,47, 0,9,190,
    0,8,15, 0,8,143, 0,8,79, 0,9,254,
    96,7,256, 0,8,80, 0,8,16, 84,8,115,
    82,7,31, 0,8,112, 0,8,48, 0,9,193,

    80,7,10, 0,8,96, 0,8,32, 0,9,161,
    0,8,0, 0,8,128, 0,8,64, 0,9,225,
    80,7,6, 0,8,88, 0,8,24, 0,9,145,
    83,7,59, 0,8,120, 0,8,56, 0,9,209,
    81,7,17, 0,8,104, 0,8,40, 0,9,177,
    0,8,8, 0,8,136, 0,8,72, 0,9,241,
    80,7,4, 0,8,84, 0,8,20, 85,8,227,
    83,7,43, 0,8,116, 0,8,52, 0,9,201,
    81,7,13, 0,8,100, 0,8,36, 0,9,169,
    0,8,4, 0,8,132, 0,8,68, 0,9,233,
    80,7,8, 0,8,92, 0,8,28, 0,9,153,
    84,7,83, 0,8,124, 0,8,60, 0,9,217,
    82,7,23, 0,8,108, 0,8,44, 0,9,185,
    0,8,12, 0,8,140, 0,8,76, 0,9,249,
    80,7,3, 0,8,82, 0,8,18, 85,8,163,
    83,7,35, 0,8,114, 0,8,50, 0,9,197,
    81,7,11, 0,8,98, 0,8,34, 0,9,165,
    0,8,2, 0,8,130, 0,8,66, 0,9,229,
    80,7,7, 0,8,90, 0,8,26, 0,9,149,
    84,7,67, 0,8,122, 0,8,58, 0,9,213,
    82,7,19, 0,8,106, 0,8,42, 0,9,181,
    0,8,10, 0,8,138, 0,8,74, 0,9,245,
    80,7,5, 0,8,86, 0,8,22, 192,8,0,
    83,7,51, 0,8,118, 0,8,54, 0,9,205,
    81,7,15, 0,8,102, 0,8,38, 0,9,173,
    0,8,6, 0,8,134, 0,8,70, 0,9,237,
    80,7,9, 0,8,94, 0,8,30, 0,9,157,
    84,7,99, 0,8,126, 0,8,62, 0,9,221,
    82,7,27, 0,8,110, 0,8,46, 0,9,189,
    0,8,14, 0,8,142, 0,8,78, 0,9,253,
    96,7,256, 0,8,81, 0,8,17, 85,8,131,
    82,7,31, 0,8,113, 0,8,49, 0,9,195,
    80,7,10, 0,8,97, 0,8,33, 0,9,163,
    0,8,1, 0,8,129, 0,8,65, 0,9,227,
    80,7,6, 0,8,89, 0,8,25, 0,9,147,
    83,7,59, 0,8,121, 0,8,57, 0,9,211,
    81,7,17, 0,8,105, 0,8,41, 0,9,179,
    0,8,9, 0,8,137, 0,8,73, 0,9,243,
    80,7,4, 0,8,85, 0,8,21, 80,8,258,
    83,7,43, 0,8,117, 0,8,53, 0,9,203,
    81,7,13, 0,8,101, 0,8,37, 0,9,171,
    0,8,5, 0,8,133, 0,8,69, 0,9,235,
    80,7,8, 0,8,93, 0,8,29, 0,9,155,
    84,7,83, 0,8,125, 0,8,61, 0,9,219,
    82,7,23, 0,8,109, 0,8,45, 0,9,187,
    0,8,13, 0,8,141, 0,8,77, 0,9,251,
    80,7,3, 0,8,83, 0,8,19, 85,8,195,
    83,7,35, 0,8,115, 0,8,51, 0,9,199,
    81,7,11, 0,8,99, 0,8,35, 0,9,167,
    0,8,3, 0,8,131, 0,8,67, 0,9,231,
    80,7,7, 0,8,91, 0,8,27, 0,9,151,
    84,7,67, 0,8,123, 0,8,59, 0,9,215,
    82,7,19, 0,8,107, 0,8,43, 0,9,183,
    0,8,11, 0,8,139, 0,8,75, 0,9,247,
    80,7,5, 0,8,87, 0,8,23, 192,8,0,
    83,7,51, 0,8,119, 0,8,55, 0,9,207,
    81,7,15, 0,8,103, 0,8,39, 0,9,175,
    0,8,7, 0,8,135, 0,8,71, 0,9,239,
    80,7,9, 0,8,95, 0,8,31, 0,9,159,
    84,7,99, 0,8,127, 0,8,63, 0,9,223,
    82,7,27, 0,8,111, 0,8,47, 0,9,191,
    0,8,15, 0,8,143, 0,8,79, 0,9,255
];
var fixed_td = [
    80,5,1, 87,5,257, 83,5,17, 91,5,4097,
    81,5,5, 89,5,1025, 85,5,65, 93,5,16385,
    80,5,3, 88,5,513, 84,5,33, 92,5,8193,
    82,5,9, 90,5,2049, 86,5,129, 192,5,24577,
    80,5,2, 87,5,385, 83,5,25, 91,5,6145,
    81,5,7, 89,5,1537, 85,5,97, 93,5,24577,
    80,5,4, 88,5,769, 84,5,49, 92,5,12289,
    82,5,13, 90,5,3073, 86,5,193, 192,5,24577
];

  // Tables for deflate from PKZIP's appnote.txt.
  var cplens = [ // Copy lengths for literal codes 257..285
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258, 0, 0
  ];

  // see note #13 above about 258
  var cplext = [ // Extra bits for literal codes 257..285
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0, 112, 112  // 112==invalid
  ];

 var cpdist = [ // Copy offsets for distance codes 0..29
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
        8193, 12289, 16385, 24577
  ];

  var cpdext = [ // Extra bits for distance codes
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11,
        12, 12, 13, 13];

//
// ZStream.java
//

function ZStream() {
}


ZStream.prototype.inflateInit = function(w, nowrap) {
    if (!w) {
	w = DEF_WBITS;
    }
    if (nowrap) {
	nowrap = false;
    }
    this.istate = new Inflate();
    return this.istate.inflateInit(this, nowrap?-w:w);
}

ZStream.prototype.inflate = function(f) {
    if(this.istate==null) return Z_STREAM_ERROR;
    return this.istate.inflate(this, f);
}

ZStream.prototype.inflateEnd = function(){
    if(this.istate==null) return Z_STREAM_ERROR;
    var ret=istate.inflateEnd(this);
    this.istate = null;
    return ret;
}
ZStream.prototype.inflateSync = function(){
    // if(istate == null) return Z_STREAM_ERROR;
    return istate.inflateSync(this);
}
ZStream.prototype.inflateSetDictionary = function(dictionary, dictLength){
    // if(istate == null) return Z_STREAM_ERROR;
    return istate.inflateSetDictionary(this, dictionary, dictLength);
}

/*

  public int deflateInit(int level){
    return deflateInit(level, MAX_WBITS);
  }
  public int deflateInit(int level, boolean nowrap){
    return deflateInit(level, MAX_WBITS, nowrap);
  }
  public int deflateInit(int level, int bits){
    return deflateInit(level, bits, false);
  }
  public int deflateInit(int level, int bits, boolean nowrap){
    dstate=new Deflate();
    return dstate.deflateInit(this, level, nowrap?-bits:bits);
  }
  public int deflate(int flush){
    if(dstate==null){
      return Z_STREAM_ERROR;
    }
    return dstate.deflate(this, flush);
  }
  public int deflateEnd(){
    if(dstate==null) return Z_STREAM_ERROR;
    int ret=dstate.deflateEnd();
    dstate=null;
    return ret;
  }
  public int deflateParams(int level, int strategy){
    if(dstate==null) return Z_STREAM_ERROR;
    return dstate.deflateParams(this, level, strategy);
  }
  public int deflateSetDictionary (byte[] dictionary, int dictLength){
    if(dstate == null)
      return Z_STREAM_ERROR;
    return dstate.deflateSetDictionary(this, dictionary, dictLength);
  }

*/

/*
  // Flush as much pending output as possible. All deflate() output goes
  // through this function so some applications may wish to modify it
  // to avoid allocating a large strm->next_out buffer and copying into it.
  // (See also read_buf()).
  void flush_pending(){
    int len=dstate.pending;

    if(len>avail_out) len=avail_out;
    if(len==0) return;

    if(dstate.pending_buf.length<=dstate.pending_out ||
       next_out.length<=next_out_index ||
       dstate.pending_buf.length<(dstate.pending_out+len) ||
       next_out.length<(next_out_index+len)){
      System.out.println(dstate.pending_buf.length+", "+dstate.pending_out+
			 ", "+next_out.length+", "+next_out_index+", "+len);
      System.out.println("avail_out="+avail_out);
    }

    System.arraycopy(dstate.pending_buf, dstate.pending_out,
		     next_out, next_out_index, len);

    next_out_index+=len;
    dstate.pending_out+=len;
    total_out+=len;
    avail_out-=len;
    dstate.pending-=len;
    if(dstate.pending==0){
      dstate.pending_out=0;
    }
  }

  // Read a new buffer from the current input stream, update the adler32
  // and total number of bytes read.  All deflate() input goes through
  // this function so some applications may wish to modify it to avoid
  // allocating a large strm->next_in buffer and copying from it.
  // (See also flush_pending()).
  int read_buf(byte[] buf, int start, int size) {
    int len=avail_in;

    if(len>size) len=size;
    if(len==0) return 0;

    avail_in-=len;

    if(dstate.noheader==0) {
      adler=_adler.adler32(adler, next_in, next_in_index, len);
    }
    System.arraycopy(next_in, next_in_index, buf, start, len);
    next_in_index  += len;
    total_in += len;
    return len;
  }

  public void free(){
    next_in=null;
    next_out=null;
    msg=null;
    _adler=null;
  }
}
*/


//
// Inflate.java
//

function Inflate() {
    this.was = [0];
}

Inflate.prototype.inflateReset = function(z) {
    if(z == null || z.istate == null) return Z_STREAM_ERROR;
    
    z.total_in = z.total_out = 0;
    z.msg = null;
    z.istate.mode = z.istate.nowrap!=0 ? BLOCKS : METHOD;
    z.istate.blocks.reset(z, null);
    return Z_OK;
}

Inflate.prototype.inflateEnd = function(z){
    if(this.blocks != null)
      this.blocks.free(z);
    this.blocks=null;
    return Z_OK;
}

Inflate.prototype.inflateInit = function(z, w){
    z.msg = null;
    this.blocks = null;

    // handle undocumented nowrap option (no zlib header or check)
    nowrap = 0;
    if(w < 0){
      w = - w;
      nowrap = 1;
    }

    // set window size
    if(w<8 ||w>15){
      this.inflateEnd(z);
      return Z_STREAM_ERROR;
    }
    this.wbits=w;

    z.istate.blocks=new InfBlocks(z, 
				  z.istate.nowrap!=0 ? null : this,
				  1<<w);

    // reset state
    this.inflateReset(z);
    return Z_OK;
  }

Inflate.prototype.inflate = function(z, f){
    var r, b;

    if(z == null || z.istate == null || z.next_in == null)
      return Z_STREAM_ERROR;
    f = f == Z_FINISH ? Z_BUF_ERROR : Z_OK;
    r = Z_BUF_ERROR;
    while (true){
      switch (z.istate.mode){
      case METHOD:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        if(((z.istate.method = z.next_in[z.next_in_index++])&0xf)!=Z_DEFLATED){
          z.istate.mode = BAD;
          z.msg="unknown compression method";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }
        if((z.istate.method>>4)+8>z.istate.wbits){
          z.istate.mode = BAD;
          z.msg="invalid window size";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }
        z.istate.mode=FLAG;
      case FLAG:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        b = (z.next_in[z.next_in_index++])&0xff;

        if((((z.istate.method << 8)+b) % 31)!=0){
          z.istate.mode = BAD;
          z.msg = "incorrect header check";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }

        if((b&PRESET_DICT)==0){
          z.istate.mode = BLOCKS;
          break;
        }
        z.istate.mode = DICT4;
      case DICT4:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need=((z.next_in[z.next_in_index++]&0xff)<<24)&0xff000000;
        z.istate.mode=DICT3;
      case DICT3:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<16)&0xff0000;
        z.istate.mode=DICT2;
      case DICT2:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<8)&0xff00;
        z.istate.mode=DICT1;
      case DICT1:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need += (z.next_in[z.next_in_index++]&0xff);
        z.adler = z.istate.need;
        z.istate.mode = DICT0;
        return Z_NEED_DICT;
      case DICT0:
        z.istate.mode = BAD;
        z.msg = "need dictionary";
        z.istate.marker = 0;       // can try inflateSync
        return Z_STREAM_ERROR;
      case BLOCKS:

        r = z.istate.blocks.proc(z, r);
        if(r == Z_DATA_ERROR){
          z.istate.mode = BAD;
          z.istate.marker = 0;     // can try inflateSync
          break;
        }
        if(r == Z_OK){
          r = f;
        }
        if(r != Z_STREAM_END){
          return r;
        }
        r = f;
        z.istate.blocks.reset(z, z.istate.was);
        if(z.istate.nowrap!=0){
          z.istate.mode=DONE;
          break;
        }
        z.istate.mode=CHECK4;
      case CHECK4:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need=((z.next_in[z.next_in_index++]&0xff)<<24)&0xff000000;
        z.istate.mode=CHECK3;
      case CHECK3:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<16)&0xff0000;
        z.istate.mode = CHECK2;
      case CHECK2:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<8)&0xff00;
        z.istate.mode = CHECK1;
      case CHECK1:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=(z.next_in[z.next_in_index++]&0xff);

        if(((z.istate.was[0])) != ((z.istate.need))){
          z.istate.mode = BAD;
          z.msg = "incorrect data check";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }

        z.istate.mode = DONE;
      case DONE:
        return Z_STREAM_END;
      case BAD:
        return Z_DATA_ERROR;
      default:
        return Z_STREAM_ERROR;
      }
    }
  }


Inflate.prototype.inflateSetDictionary = function(z,  dictionary, dictLength) {
    var index=0;
    var length = dictLength;
    if(z==null || z.istate == null|| z.istate.mode != DICT0)
      return Z_STREAM_ERROR;

    if(z._adler.adler32(1, dictionary, 0, dictLength)!=z.adler){
      return Z_DATA_ERROR;
    }

    z.adler = z._adler.adler32(0, null, 0, 0);

    if(length >= (1<<z.istate.wbits)){
      length = (1<<z.istate.wbits)-1;
      index=dictLength - length;
    }
    z.istate.blocks.set_dictionary(dictionary, index, length);
    z.istate.mode = BLOCKS;
    return Z_OK;
  }

//  static private byte[] mark = {(byte)0, (byte)0, (byte)0xff, (byte)0xff};
var mark = [0, 0, 255, 255]

Inflate.prototype.inflateSync = function(z){
    var n;       // number of bytes to look at
    var p;       // pointer to bytes
    var m;       // number of marker bytes found in a row
    var r, w;   // temporaries to save total_in and total_out

    // set up
    if(z == null || z.istate == null)
      return Z_STREAM_ERROR;
    if(z.istate.mode != BAD){
      z.istate.mode = BAD;
      z.istate.marker = 0;
    }
    if((n=z.avail_in)==0)
      return Z_BUF_ERROR;
    p=z.next_in_index;
    m=z.istate.marker;

    // search
    while (n!=0 && m < 4){
      if(z.next_in[p] == mark[m]){
        m++;
      }
      else if(z.next_in[p]!=0){
        m = 0;
      }
      else{
        m = 4 - m;
      }
      p++; n--;
    }

    // restore
    z.total_in += p-z.next_in_index;
    z.next_in_index = p;
    z.avail_in = n;
    z.istate.marker = m;

    // return no joy or set up to restart on a new block
    if(m != 4){
      return Z_DATA_ERROR;
    }
    r=z.total_in;  w=z.total_out;
    this.inflateReset(z);
    z.total_in=r;  z.total_out = w;
    z.istate.mode = BLOCKS;
    return Z_OK;
}

  // Returns true if inflate is currently at the end of a block generated
  // by Z_SYNC_FLUSH or Z_FULL_FLUSH. This function is used by one PPP
  // implementation to provide an additional safety check. PPP uses Z_SYNC_FLUSH
  // but removes the length bytes of the resulting empty stored block. When
  // decompressing, PPP checks that at the end of input packet, inflate is
  // waiting for these length bytes.
Inflate.prototype.inflateSyncPoint = function(z){
    if(z == null || z.istate == null || z.istate.blocks == null)
      return Z_STREAM_ERROR;
    return z.istate.blocks.sync_point();
}


//
// InfBlocks.java
//

var INFBLOCKS_BORDER = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

function InfBlocks(z, checkfn, w) {
    this.hufts=new Int32Array(MANY*3);
    this.window=new Uint8Array(w);
    this.end=w;
    this.checkfn = checkfn;
    this.mode = IB_TYPE;
    this.reset(z, null);

    this.left = 0;            // if STORED, bytes left to copy 

    this.table = 0;           // table lengths (14 bits) 
    this.index = 0;           // index into blens (or border) 
    this.blens = null;         // bit lengths of codes 
    this.bb=new Int32Array(1); // bit length tree depth 
    this.tb=new Int32Array(1); // bit length decoding tree 

    this.codes = new InfCodes();

    this.last = 0;            // true if this block is the last block 

  // mode independent information 
    this.bitk = 0;            // bits in bit buffer 
    this.bitb = 0;            // bit buffer 
    this.read = 0;            // window read pointer 
    this.write = 0;           // window write pointer 
    this.check = 0;          // check on output 

    this.inftree=new InfTree();
}




InfBlocks.prototype.reset = function(z, c){
    if(c) c[0]=this.check;
    if(this.mode==IB_CODES){
      this.codes.free(z);
    }
    this.mode=IB_TYPE;
    this.bitk=0;
    this.bitb=0;
    this.read=this.write=0;

    if(this.checkfn)
      z.adler=this.check=z._adler.adler32(0, null, 0, 0);
  }

 InfBlocks.prototype.proc = function(z, r){
    var t;              // temporary storage
    var b;              // bit buffer
    var k;              // bits in bit buffer
    var p;              // input data pointer
    var n;              // bytes available there
    var q;              // output window write pointer
    var m;              // bytes to end of window or read pointer

    // copy input/output information to locals (UPDATE macro restores)
    {p=z.next_in_index;n=z.avail_in;b=this.bitb;k=this.bitk;}
    {q=this.write;m=(q<this.read ? this.read-q-1 : this.end-q);}

    // process input based on current state
    while(true){
      switch (this.mode){
      case IB_TYPE:

	while(k<(3)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}
	t = (b & 7);
	this.last = t & 1;

	switch (t >>> 1){
        case 0:                         // stored 
          {b>>>=(3);k-=(3);}
          t = k & 7;                    // go to byte boundary

          {b>>>=(t);k-=(t);}
          this.mode = IB_LENS;                  // get length of stored block
          break;
        case 1:                         // fixed
          {
              var bl=new Int32Array(1);
	      var bd=new Int32Array(1);
              var tl=[];
	      var td=[];

	      inflate_trees_fixed(bl, bd, tl, td, z);
              this.codes.init(bl[0], bd[0], tl[0], 0, td[0], 0, z);
          }

          {b>>>=(3);k-=(3);}

          this.mode = IB_CODES;
          break;
        case 2:                         // dynamic

          {b>>>=(3);k-=(3);}

          this.mode = IB_TABLE;
          break;
        case 3:                         // illegal

          {b>>>=(3);k-=(3);}
          this.mode = BAD;
          z.msg = "invalid block type";
          r = Z_DATA_ERROR;

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z,r);
	}
	break;
      case IB_LENS:
	while(k<(32)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	if ((((~b) >>> 16) & 0xffff) != (b & 0xffff)){
	  this.mode = BAD;
	  z.msg = "invalid stored block lengths";
	  r = Z_DATA_ERROR;

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z,r);
	}
	this.left = (b & 0xffff);
	b = k = 0;                       // dump bits
	this.mode = this.left!=0 ? IB_STORED : (this.last!=0 ? IB_DRY : IB_TYPE);
	break;
      case IB_STORED:
	if (n == 0){
	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  write=q;
	  return this.inflate_flush(z,r);
	}

	if(m==0){
	  if(q==end&&read!=0){
	    q=0; m=(q<this.read ? this.read-q-1 : this.end-q);
	  }
	  if(m==0){
	    this.write=q; 
	    r=this.inflate_flush(z,r);
	    q=this.write; m = (q < this.read ? this.read-q-1 : this.end-q);
	    if(q==this.end && this.read != 0){
	      q=0; m = (q < this.read ? this.read-q-1 : this.end-q);
	    }
	    if(m==0){
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    }
	  }
	}
	r=Z_OK;

	t = this.left;
	if(t>n) t = n;
	if(t>m) t = m;
	arrayCopy(z.next_in, p, window, q, t);
	p += t;  n -= t;
	q += t;  m -= t;
	if ((this.left -= t) != 0)
	  break;
	this.mode = (this.last != 0 ? IB_DRY : IB_TYPE);
	break;
      case IB_TABLE:

	while(k<(14)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.table = t = (b & 0x3fff);
	if ((t & 0x1f) > 29 || ((t >> 5) & 0x1f) > 29)
	  {
	    this.mode = IB_BAD;
	    z.msg = "too many length or distance symbols";
	    r = Z_DATA_ERROR;

	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  }
	t = 258 + (t & 0x1f) + ((t >> 5) & 0x1f);
	if(this.blens==null || this.blens.length<t){
	    this.blens=new Int32Array(t);
	}
	else{
	  for(var i=0; i<t; i++){
              this.blens[i]=0;
          }
	}

	{b>>>=(14);k-=(14);}

	this.index = 0;
	mode = IB_BTREE;
      case IB_BTREE:
	while (this.index < 4 + (this.table >>> 10)){
	  while(k<(3)){
	    if(n!=0){
	      r=Z_OK;
	    }
	    else{
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;
	      z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    };
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;
	    k+=8;
	  }

	  this.blens[INFBLOCKS_BORDER[this.index++]] = b&7;

	  {b>>>=(3);k-=(3);}
	}

	while(this.index < 19){
	  this.blens[INFBLOCKS_BORDER[this.index++]] = 0;
	}

	this.bb[0] = 7;
	t = this.inftree.inflate_trees_bits(this.blens, this.bb, this.tb, this.hufts, z);
	if (t != Z_OK){
	  r = t;
	  if (r == Z_DATA_ERROR){
	    this.blens=null;
	    this.mode = IB_BAD;
	  }

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  write=q;
	  return this.inflate_flush(z,r);
	}

	this.index = 0;
	this.mode = IB_DTREE;
      case IB_DTREE:
	while (true){
	  t = this.table;
	  if(!(this.index < 258 + (t & 0x1f) + ((t >> 5) & 0x1f))){
	    break;
	  }

	  var h; //int[]
	  var i, j, c;

	  t = this.bb[0];

	  while(k<(t)){
	    if(n!=0){
	      r=Z_OK;
	    }
	    else{
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;
	      z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    };
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;
	    k+=8;
	  }

//	  if (this.tb[0]==-1){
//            dlog("null...");
//	  }

	  t=this.hufts[(this.tb[0]+(b & inflate_mask[t]))*3+1];
	  c=this.hufts[(this.tb[0]+(b & inflate_mask[t]))*3+2];

	  if (c < 16){
	    b>>>=(t);k-=(t);
	    this.blens[this.index++] = c;
	  }
	  else { // c == 16..18
	    i = c == 18 ? 7 : c - 14;
	    j = c == 18 ? 11 : 3;

	    while(k<(t+i)){
	      if(n!=0){
		r=Z_OK;
	      }
	      else{
		this.bitb=b; this.bitk=k; 
		z.avail_in=n;
		z.total_in+=p-z.next_in_index;z.next_in_index=p;
		this.write=q;
		return this.inflate_flush(z,r);
	      };
	      n--;
	      b|=(z.next_in[p++]&0xff)<<k;
	      k+=8;
	    }

	    b>>>=(t);k-=(t);

	    j += (b & inflate_mask[i]);

	    b>>>=(i);k-=(i);

	    i = this.index;
	    t = this.table;
	    if (i + j > 258 + (t & 0x1f) + ((t >> 5) & 0x1f) ||
		(c == 16 && i < 1)){
	      this.blens=null;
	      this.mode = IB_BAD;
	      z.msg = "invalid bit length repeat";
	      r = Z_DATA_ERROR;

	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    }

	    c = c == 16 ? this.blens[i-1] : 0;
	    do{
	      this.blens[i++] = c;
	    }
	    while (--j!=0);
	    this.index = i;
	  }
	}

	this.tb[0]=-1;
	{
	    var bl=new Int32Array(1);
	    var bd=new Int32Array(1);
	    var tl=new Int32Array(1);
	    var td=new Int32Array(1);
	    bl[0] = 9;         // must be <= 9 for lookahead assumptions
	    bd[0] = 6;         // must be <= 9 for lookahead assumptions

	    t = this.table;
	    t = this.inftree.inflate_trees_dynamic(257 + (t & 0x1f), 
					      1 + ((t >> 5) & 0x1f),
					      this.blens, bl, bd, tl, td, this.hufts, z);

	    if (t != Z_OK){
	        if (t == Z_DATA_ERROR){
	            this.blens=null;
	            this.mode = BAD;
	        }
	        r = t;

	        this.bitb=b; this.bitk=k; 
	        z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	        this.write=q;
	        return this.inflate_flush(z,r);
	    }
	    this.codes.init(bl[0], bd[0], this.hufts, tl[0], this.hufts, td[0], z);
	}
	this.mode = IB_CODES;
      case IB_CODES:
	this.bitb=b; this.bitk=k;
	z.avail_in=n; z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;

	if ((r = this.codes.proc(this, z, r)) != Z_STREAM_END){
	  return this.inflate_flush(z, r);
	}
	r = Z_OK;
	this.codes.free(z);

	p=z.next_in_index; n=z.avail_in;b=this.bitb;k=this.bitk;
	q=this.write;m = (q < this.read ? this.read-q-1 : this.end-q);

	if (this.last==0){
	  this.mode = IB_TYPE;
	  break;
	}
	this.mode = IB_DRY;
      case IB_DRY:
	this.write=q; 
	r = this.inflate_flush(z, r); 
	q=this.write; m = (q < this.read ? this.read-q-1 : this.end-q);
	if (this.read != this.write){
	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z, r);
	}
	mode = DONE;
      case IB_DONE:
	r = Z_STREAM_END;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);
      case IB_BAD:
	r = Z_DATA_ERROR;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);

      default:
	r = Z_STREAM_ERROR;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);
      }
    }
  }

InfBlocks.prototype.free = function(z){
    this.reset(z, null);
    this.window=null;
    this.hufts=null;
}

InfBlocks.prototype.set_dictionary = function(d, start, n){
    arrayCopy(d, start, window, 0, n);
    this.read = this.write = n;
}

  // Returns true if inflate is currently at the end of a block generated
  // by Z_SYNC_FLUSH or Z_FULL_FLUSH. 
InfBlocks.prototype.sync_point = function(){
    return this.mode == IB_LENS;
}

  // copy as much as possible from the sliding window to the output area
InfBlocks.prototype.inflate_flush = function(z, r){
    var n;
    var p;
    var q;

    // local copies of source and destination pointers
    p = z.next_out_index;
    q = this.read;

    // compute number of bytes to copy as far as end of window
    n = ((q <= this.write ? this.write : this.end) - q);
    if (n > z.avail_out) n = z.avail_out;
    if (n!=0 && r == Z_BUF_ERROR) r = Z_OK;

    // update counters
    z.avail_out -= n;
    z.total_out += n;

    // update check information
    if(this.checkfn != null)
      z.adler=this.check=z._adler.adler32(this.check, this.window, q, n);

    // copy as far as end of window
    arrayCopy(this.window, q, z.next_out, p, n);
    p += n;
    q += n;

    // see if more to copy at beginning of window
    if (q == this.end){
      // wrap pointers
      q = 0;
      if (this.write == this.end)
        this.write = 0;

      // compute bytes to copy
      n = this.write - q;
      if (n > z.avail_out) n = z.avail_out;
      if (n!=0 && r == Z_BUF_ERROR) r = Z_OK;

      // update counters
      z.avail_out -= n;
      z.total_out += n;

      // update check information
      if(this.checkfn != null)
	z.adler=this.check=z._adler.adler32(this.check, this.window, q, n);

      // copy
      arrayCopy(this.window, q, z.next_out, p, n);
      p += n;
      q += n;
    }

    // update pointers
    z.next_out_index = p;
    this.read = q;

    // done
    return r;
  }

//
// InfCodes.java
//

var IC_START=0;  // x: set up for LEN
var IC_LEN=1;    // i: get length/literal/eob next
var IC_LENEXT=2; // i: getting length extra (have base)
var IC_DIST=3;   // i: get distance next
var IC_DISTEXT=4;// i: getting distance extra
var IC_COPY=5;   // o: copying bytes in window, waiting for space
var IC_LIT=6;    // o: got literal, waiting for output space
var IC_WASH=7;   // o: got eob, possibly still output waiting
var IC_END=8;    // x: got eob and all data flushed
var IC_BADCODE=9;// x: got error

function InfCodes() {
}

InfCodes.prototype.init = function(bl, bd, tl, tl_index, td, td_index, z) {
    this.mode=IC_START;
    this.lbits=bl;
    this.dbits=bd;
    this.ltree=tl;
    this.ltree_index=tl_index;
    this.dtree = td;
    this.dtree_index=td_index;
    this.tree=null;
}

InfCodes.prototype.proc = function(s, z, r){ 
    var j;              // temporary storage
    var t;              // temporary pointer (int[])
    var tindex;         // temporary pointer
    var e;              // extra bits or operation
    var b=0;            // bit buffer
    var k=0;            // bits in bit buffer
    var p=0;            // input data pointer
    var n;              // bytes available there
    var q;              // output window write pointer
    var m;              // bytes to end of window or read pointer
    var f;              // pointer to copy strings from

    // copy input/output information to locals (UPDATE macro restores)
    p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

    // process input and output based on current state
    while (true){
      switch (this.mode){
	// waiting for "i:"=input, "o:"=output, "x:"=nothing
      case IC_START:         // x: set up for LEN
	if (m >= 258 && n >= 10){

	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;
	  r = this.inflate_fast(this.lbits, this.dbits, 
			   this.ltree, this.ltree_index, 
			   this.dtree, this.dtree_index,
			   s, z);

	  p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
	  q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	  if (r != Z_OK){
	    this.mode = r == Z_STREAM_END ? IC_WASH : IC_BADCODE;
	    break;
	  }
	}
	this.need = this.lbits;
	this.tree = this.ltree;
	this.tree_index=this.ltree_index;

	this.mode = IC_LEN;
      case IC_LEN:           // i: get length/literal/eob next
	j = this.need;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	tindex=(this.tree_index+(b&inflate_mask[j]))*3;

	b>>>=(this.tree[tindex+1]);
	k-=(this.tree[tindex+1]);

	e=this.tree[tindex];

	if(e == 0){               // literal
	  this.lit = this.tree[tindex+2];
	  this.mode = IC_LIT;
	  break;
	}
	if((e & 16)!=0 ){          // length
	  this.get = e & 15;
	  this.len = this.tree[tindex+2];
	  this.mode = IC_LENEXT;
	  break;
	}
	if ((e & 64) == 0){        // next table
	  this.need = e;
	  this.tree_index = tindex/3 + this.tree[tindex+2];
	  break;
	}
	if ((e & 32)!=0){               // end of block
	  this.mode = IC_WASH;
	  break;
	}
	this.mode = IC_BADCODE;        // invalid code
	z.msg = "invalid literal/length code";
	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_LENEXT:        // i: getting length extra (have base)
	j = this.get;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.len += (b & inflate_mask[j]);

	b>>=j;
	k-=j;

	this.need = this.dbits;
	this.tree = this.dtree;
	this.tree_index = this.dtree_index;
	this.mode = IC_DIST;
      case IC_DIST:          // i: get distance next
	j = this.need;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	tindex=(this.tree_index+(b & inflate_mask[j]))*3;

	b>>=this.tree[tindex+1];
	k-=this.tree[tindex+1];

	e = (this.tree[tindex]);
	if((e & 16)!=0){               // distance
	  this.get = e & 15;
	  this.dist = this.tree[tindex+2];
	  this.mode = IC_DISTEXT;
	  break;
	}
	if ((e & 64) == 0){        // next table
	  this.need = e;
	  this.tree_index = tindex/3 + this.tree[tindex+2];
	  break;
	}
	this.mode = IC_BADCODE;        // invalid code
	z.msg = "invalid distance code";
	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_DISTEXT:       // i: getting distance extra
	j = this.get;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.dist += (b & inflate_mask[j]);

	b>>=j;
	k-=j;

	this.mode = IC_COPY;
      case IC_COPY:          // o: copying bytes in window, waiting for space
        f = q - this.dist;
        while(f < 0){     // modulo window size-"while" instead
          f += s.end;     // of "if" handles invalid distances
	}
	while (this.len!=0){

	  if(m==0){
	    if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	    if(m==0){
	      s.write=q; r=s.inflate_flush(z,r);
	      q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	      if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}

	      if(m==0){
		s.bitb=b;s.bitk=k;
		z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
		s.write=q;
		return s.inflate_flush(z,r);
	      }  
	    }
	  }

	  s.window[q++]=s.window[f++]; m--;

	  if (f == s.end)
            f = 0;
	  this.len--;
	}
	this.mode = IC_START;
	break;
      case IC_LIT:           // o: got literal, waiting for output space
	if(m==0){
	  if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	  if(m==0){
	    s.write=q; r=s.inflate_flush(z,r);
	    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	    if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	    if(m==0){
	      s.bitb=b;s.bitk=k;
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      s.write=q;
	      return s.inflate_flush(z,r);
	    }
	  }
	}
	r=Z_OK;

	s.window[q++]=this.lit; m--;

	this.mode = IC_START;
	break;
      case IC_WASH:           // o: got eob, possibly more output
	if (k > 7){        // return unused byte, if any
	  k -= 8;
	  n++;
	  p--;             // can always return one
	}

	s.write=q; r=s.inflate_flush(z,r);
	q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	if (s.read != s.write){
	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;
	  return s.inflate_flush(z,r);
	}
	this.mode = IC_END;
      case IC_END:
	r = Z_STREAM_END;
	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_BADCODE:       // x: got error

	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      default:
	r = Z_STREAM_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);
      }
    }
  }

InfCodes.prototype.free = function(z){
    //  ZFREE(z, c);
}

  // Called with number of bytes left to write in window at least 258
  // (the maximum string length) and number of input bytes available
  // at least ten.  The ten bytes are six bytes for the longest length/
  // distance pair plus four bytes for overloading the bit buffer.

InfCodes.prototype.inflate_fast = function(bl, bd, tl, tl_index, td, td_index, s, z) {
    var t;                // temporary pointer
    var   tp;             // temporary pointer (int[])
    var tp_index;         // temporary pointer
    var e;                // extra bits or operation
    var b;                // bit buffer
    var k;                // bits in bit buffer
    var p;                // input data pointer
    var n;                // bytes available there
    var q;                // output window write pointer
    var m;                // bytes to end of window or read pointer
    var ml;               // mask for literal/length tree
    var md;               // mask for distance tree
    var c;                // bytes to copy
    var d;                // distance back to copy from
    var r;                // copy source pointer

    var tp_index_t_3;     // (tp_index+t)*3

    // load input, output, bit values
    p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

    // initialize masks
    ml = inflate_mask[bl];
    md = inflate_mask[bd];

    // do until not enough input or output space for fast loop
    do {                          // assume called with m >= 258 && n >= 10
      // get literal/length code
      while(k<(20)){              // max bits for literal/length code
	n--;
	b|=(z.next_in[p++]&0xff)<<k;k+=8;
      }

      t= b&ml;
      tp=tl; 
      tp_index=tl_index;
      tp_index_t_3=(tp_index+t)*3;
      if ((e = tp[tp_index_t_3]) == 0){
	b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	s.window[q++] = tp[tp_index_t_3+2];
	m--;
	continue;
      }
      do {

	b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	if((e&16)!=0){
	  e &= 15;
	  c = tp[tp_index_t_3+2] + (b & inflate_mask[e]);

	  b>>=e; k-=e;

	  // decode distance base of block to copy
	  while(k<(15)){           // max bits for distance code
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;k+=8;
	  }

	  t= b&md;
	  tp=td;
	  tp_index=td_index;
          tp_index_t_3=(tp_index+t)*3;
	  e = tp[tp_index_t_3];

	  do {

	    b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	    if((e&16)!=0){
	      // get extra bits to add to distance base
	      e &= 15;
	      while(k<(e)){         // get extra bits (up to 13)
		n--;
		b|=(z.next_in[p++]&0xff)<<k;k+=8;
	      }

	      d = tp[tp_index_t_3+2] + (b&inflate_mask[e]);

	      b>>=(e); k-=(e);

	      // do the copy
	      m -= c;
	      if (q >= d){                // offset before dest
		//  just copy
		r=q-d;
		if(q-r>0 && 2>(q-r)){           
		  s.window[q++]=s.window[r++]; // minimum count is three,
		  s.window[q++]=s.window[r++]; // so unroll loop a little
		  c-=2;
		}
		else{
		  s.window[q++]=s.window[r++]; // minimum count is three,
		  s.window[q++]=s.window[r++]; // so unroll loop a little
		  c-=2;
		}
	      }
	      else{                  // else offset after destination
                r=q-d;
                do{
                  r+=s.end;          // force pointer in window
                }while(r<0);         // covers invalid distances
		e=s.end-r;
		if(c>e){             // if source crosses,
		  c-=e;              // wrapped copy
		  if(q-r>0 && e>(q-r)){           
		    do{s.window[q++] = s.window[r++];}
		    while(--e!=0);
		  }
		  else{
		    arrayCopy(s.window, r, s.window, q, e);
		    q+=e; r+=e; e=0;
		  }
		  r = 0;                  // copy rest from start of window
		}

	      }

	      // copy all or what's left
              do{s.window[q++] = s.window[r++];}
		while(--c!=0);
	      break;
	    }
	    else if((e&64)==0){
	      t+=tp[tp_index_t_3+2];
	      t+=(b&inflate_mask[e]);
	      tp_index_t_3=(tp_index+t)*3;
	      e=tp[tp_index_t_3];
	    }
	    else{
	      z.msg = "invalid distance code";

	      c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

	      s.bitb=b;s.bitk=k;
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      s.write=q;

	      return Z_DATA_ERROR;
	    }
	  }
	  while(true);
	  break;
	}

	if((e&64)==0){
	  t+=tp[tp_index_t_3+2];
	  t+=(b&inflate_mask[e]);
	  tp_index_t_3=(tp_index+t)*3;
	  if((e=tp[tp_index_t_3])==0){

	    b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	    s.window[q++]=tp[tp_index_t_3+2];
	    m--;
	    break;
	  }
	}
	else if((e&32)!=0){

	  c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;
 
	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;

	  return Z_STREAM_END;
	}
	else{
	  z.msg="invalid literal/length code";

	  c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;

	  return Z_DATA_ERROR;
	}
      } 
      while(true);
    } 
    while(m>=258 && n>= 10);

    // not enough input or output--restore pointers and return
    c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

    s.bitb=b;s.bitk=k;
    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
    s.write=q;

    return Z_OK;
}

//
// InfTree.java
//

function InfTree() {
}

InfTree.prototype.huft_build = function(b, bindex, n, s, d, e, t, m, hp, hn, v) {

    // Given a list of code lengths and a maximum table size, make a set of
    // tables to decode that set of codes.  Return Z_OK on success, Z_BUF_ERROR
    // if the given code set is incomplete (the tables are still built in this
    // case), Z_DATA_ERROR if the input is invalid (an over-subscribed set of
    // lengths), or Z_MEM_ERROR if not enough memory.

    var a;                       // counter for codes of length k
    var f;                       // i repeats in table every f entries
    var g;                       // maximum code length
    var h;                       // table level
    var i;                       // counter, current code
    var j;                       // counter
    var k;                       // number of bits in current code
    var l;                       // bits per table (returned in m)
    var mask;                    // (1 << w) - 1, to avoid cc -O bug on HP
    var p;                       // pointer into c[], b[], or v[]
    var q;                       // points to current table
    var w;                       // bits before this table == (l * h)
    var xp;                      // pointer into x
    var y;                       // number of dummy codes added
    var z;                       // number of entries in current table

    // Generate counts for each bit length

    p = 0; i = n;
    do {
      this.c[b[bindex+p]]++; p++; i--;   // assume all entries <= BMAX
    }while(i!=0);

    if(this.c[0] == n){                // null input--all zero length codes
      t[0] = -1;
      m[0] = 0;
      return Z_OK;
    }

    // Find minimum and maximum length, bound *m by those
    l = m[0];
    for (j = 1; j <= BMAX; j++)
      if(this.c[j]!=0) break;
    k = j;                        // minimum code length
    if(l < j){
      l = j;
    }
    for (i = BMAX; i!=0; i--){
      if(this.c[i]!=0) break;
    }
    g = i;                        // maximum code length
    if(l > i){
      l = i;
    }
    m[0] = l;

    // Adjust last length count to fill out codes, if needed
    for (y = 1 << j; j < i; j++, y <<= 1){
      if ((y -= this.c[j]) < 0){
        return Z_DATA_ERROR;
      }
    }
    if ((y -= this.c[i]) < 0){
      return Z_DATA_ERROR;
    }
    this.c[i] += y;

    // Generate starting offsets into the value table for each length
    this.x[1] = j = 0;
    p = 1;  xp = 2;
    while (--i!=0) {                 // note that i == g from above
      this.x[xp] = (j += this.c[p]);
      xp++;
      p++;
    }

    // Make a table of values in order of bit lengths
    i = 0; p = 0;
    do {
      if ((j = b[bindex+p]) != 0){
        this.v[this.x[j]++] = i;
      }
      p++;
    }
    while (++i < n);
    n = this.x[g];                     // set n to length of v

    // Generate the Huffman codes and for each, make the table entries
    this.x[0] = i = 0;                 // first Huffman code is zero
    p = 0;                        // grab values in bit order
    h = -1;                       // no tables yet--level -1
    w = -l;                       // bits decoded == (l * h)
    this.u[0] = 0;                     // just to keep compilers happy
    q = 0;                        // ditto
    z = 0;                        // ditto

    // go through the bit lengths (k already is bits in shortest code)
    for (; k <= g; k++){
      a = this.c[k];
      while (a--!=0){
	// here i is the Huffman code of length k bits for value *p
	// make tables up to required level
        while (k > w + l){
          h++;
          w += l;                 // previous table always l bits
	  // compute minimum size table less than or equal to l bits
          z = g - w;
          z = (z > l) ? l : z;        // table size upper limit
          if((f=1<<(j=k-w))>a+1){     // try a k-w bit table
                                      // too few codes for k-w bit table
            f -= a + 1;               // deduct codes from patterns left
            xp = k;
            if(j < z){
              while (++j < z){        // try smaller tables up to z bits
                if((f <<= 1) <= this.c[++xp])
                  break;              // enough codes to use up j bits
                f -= this.c[xp];           // else deduct codes from patterns
              }
	    }
          }
          z = 1 << j;                 // table entries for j-bit table

	  // allocate new table
          if (this.hn[0] + z > MANY){       // (note: doesn't matter for fixed)
            return Z_DATA_ERROR;       // overflow of MANY
          }
          this.u[h] = q = /*hp+*/ this.hn[0];   // DEBUG
          this.hn[0] += z;
 
	  // connect to last table, if there is one
	  if(h!=0){
            this.x[h]=i;           // save pattern for backing up
            this.r[0]=j;     // bits in this table
            this.r[1]=l;     // bits to dump before this table
            j=i>>>(w - l);
            this.r[2] = (q - this.u[h-1] - j);               // offset to this table
            arrayCopy(this.r, 0, hp, (this.u[h-1]+j)*3, 3); // connect to last table
          }
          else{
            t[0] = q;               // first table is returned result
	  }
        }

	// set up table entry in r
        this.r[1] = (k - w);
        if (p >= n){
          this.r[0] = 128 + 64;      // out of values--invalid code
	}
        else if (v[p] < s){
          this.r[0] = (this.v[p] < 256 ? 0 : 32 + 64);  // 256 is end-of-block
          this.r[2] = this.v[p++];          // simple code is just the value
        }
        else{
          this.r[0]=(e[this.v[p]-s]+16+64); // non-simple--look up in lists
          this.r[2]=d[this.v[p++] - s];
        }

        // fill code-like entries with r
        f=1<<(k-w);
        for (j=i>>>w;j<z;j+=f){
          arrayCopy(this.r, 0, hp, (q+j)*3, 3);
	}

	// backwards increment the k-bit code i
        for (j = 1 << (k - 1); (i & j)!=0; j >>>= 1){
          i ^= j;
	}
        i ^= j;

	// backup over finished tables
        mask = (1 << w) - 1;      // needed on HP, cc -O bug
        while ((i & mask) != this.x[h]){
          h--;                    // don't need to update q
          w -= l;
          mask = (1 << w) - 1;
        }
      }
    }
    // Return Z_BUF_ERROR if we were given an incomplete table
    return y != 0 && g != 1 ? Z_BUF_ERROR : Z_OK;
}

InfTree.prototype.inflate_trees_bits = function(c, bb, tb, hp, z) {
    var result;
    this.initWorkArea(19);
    this.hn[0]=0;
    result = this.huft_build(c, 0, 19, 19, null, null, tb, bb, hp, this.hn, this.v);

    if(result == Z_DATA_ERROR){
      z.msg = "oversubscribed dynamic bit lengths tree";
    }
    else if(result == Z_BUF_ERROR || bb[0] == 0){
      z.msg = "incomplete dynamic bit lengths tree";
      result = Z_DATA_ERROR;
    }
    return result;
}

InfTree.prototype.inflate_trees_dynamic = function(nl, nd, c, bl, bd, tl, td, hp, z) {
    var result;

    // build literal/length tree
    this.initWorkArea(288);
    this.hn[0]=0;
    result = this.huft_build(c, 0, nl, 257, cplens, cplext, tl, bl, hp, this.hn, this.v);
    if (result != Z_OK || bl[0] == 0){
      if(result == Z_DATA_ERROR){
        z.msg = "oversubscribed literal/length tree";
      }
      else if (result != Z_MEM_ERROR){
        z.msg = "incomplete literal/length tree";
        result = Z_DATA_ERROR;
      }
      return result;
    }

    // build distance tree
    this.initWorkArea(288);
    result = this.huft_build(c, nl, nd, 0, cpdist, cpdext, td, bd, hp, this.hn, this.v);

    if (result != Z_OK || (bd[0] == 0 && nl > 257)){
      if (result == Z_DATA_ERROR){
        z.msg = "oversubscribed distance tree";
      }
      else if (result == Z_BUF_ERROR) {
        z.msg = "incomplete distance tree";
        result = Z_DATA_ERROR;
      }
      else if (result != Z_MEM_ERROR){
        z.msg = "empty distance tree with lengths";
        result = Z_DATA_ERROR;
      }
      return result;
    }

    return Z_OK;
}
/*
  static int inflate_trees_fixed(int[] bl,  //literal desired/actual bit depth
                                 int[] bd,  //distance desired/actual bit depth
                                 int[][] tl,//literal/length tree result
                                 int[][] td,//distance tree result 
                                 ZStream z  //for memory allocation
				 ){

*/

function inflate_trees_fixed(bl, bd, tl, td, z) {
    bl[0]=fixed_bl;
    bd[0]=fixed_bd;
    tl[0]=fixed_tl;
    td[0]=fixed_td;
    return Z_OK;
}

InfTree.prototype.initWorkArea = function(vsize){
    if(this.hn==null){
        this.hn=new Int32Array(1);
        this.v=new Int32Array(vsize);
        this.c=new Int32Array(BMAX+1);
        this.r=new Int32Array(3);
        this.u=new Int32Array(BMAX);
        this.x=new Int32Array(BMAX+1);
    }
    if(this.v.length<vsize){ 
        this.v=new Int32Array(vsize); 
    }
    for(var i=0; i<vsize; i++){this.v[i]=0;}
    for(var i=0; i<BMAX+1; i++){this.c[i]=0;}
    for(var i=0; i<3; i++){this.r[i]=0;}
//  for(int i=0; i<BMAX; i++){u[i]=0;}
    arrayCopy(this.c, 0, this.u, 0, BMAX);
//  for(int i=0; i<BMAX+1; i++){x[i]=0;}
    arrayCopy(this.c, 0, this.x, 0, BMAX+1);
}

var testArray = new Uint8Array(1);
var hasSubarray = (typeof testArray.subarray === 'function');
var hasSlice = false; /* (typeof testArray.slice === 'function'); */ // Chrome slice performance is so dire that we're currently not using it...

function arrayCopy(src, srcOffset, dest, destOffset, count) {
    if (count == 0) {
        return;
    } 
    if (!src) {
        throw "Undef src";
    } else if (!dest) {
        throw "Undef dest";
    }

    if (srcOffset == 0 && count == src.length) {
        arrayCopy_fast(src, dest, destOffset);
    } else if (hasSubarray) {
        arrayCopy_fast(src.subarray(srcOffset, srcOffset + count), dest, destOffset); 
    } else if (src.BYTES_PER_ELEMENT == 1 && count > 100) {
        arrayCopy_fast(new Uint8Array(src.buffer, src.byteOffset + srcOffset, count), dest, destOffset);
    } else { 
        arrayCopy_slow(src, srcOffset, dest, destOffset, count);
    }

}

function arrayCopy_slow(src, srcOffset, dest, destOffset, count) {

    // dlog('_slow call: srcOffset=' + srcOffset + '; destOffset=' + destOffset + '; count=' + count);

     for (var i = 0; i < count; ++i) {
        dest[destOffset + i] = src[srcOffset + i];
    }
}

function arrayCopy_fast(src, dest, destOffset) {
    dest.set(src, destOffset);
}


  // largest prime smaller than 65536
var ADLER_BASE=65521; 
  // NMAX is the largest n such that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1
var ADLER_NMAX=5552;

function adler32(adler, /* byte[] */ buf,  index, len){
    if(buf == null){ return 1; }

    var s1=adler&0xffff;
    var s2=(adler>>16)&0xffff;
    var k;

    while(len > 0) {
      k=len<ADLER_NMAX?len:ADLER_NMAX;
      len-=k;
      while(k>=16){
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        k-=16;
      }
      if(k!=0){
        do{
          s1+=buf[index++]&0xff; s2+=s1;
        }
        while(--k!=0);
      }
      s1%=ADLER_BASE;
      s2%=ADLER_BASE;
    }
    return (s2<<16)|s1;
}



function jszlib_inflate_buffer(buffer, start, length, afterUncOffset) {
    if (!start) {
        buffer = new Uint8Array(buffer);
    } else if (!length) {
        buffer = new Uint8Array(buffer, start, buffer.byteLength - start);
    } else {
        buffer = new Uint8Array(buffer, start, length);
    }

    var z = new ZStream();
    z.inflateInit(DEF_WBITS, true);
    z.next_in = buffer;
    z.next_in_index = 0;
    z.avail_in = buffer.length;

    var oBlockList = [];
    var totalSize = 0;
    while (true) {
        var obuf = new Uint8Array(32000);
        z.next_out = obuf;
        z.next_out_index = 0;
        z.avail_out = obuf.length;
        var status = z.inflate(Z_NO_FLUSH);
        if (status != Z_OK && status != Z_STREAM_END && status != Z_BUF_ERROR) {
            throw z.msg;
        }
        if (z.avail_out != 0) {
            var newob = new Uint8Array(obuf.length - z.avail_out);
            arrayCopy(obuf, 0, newob, 0, (obuf.length - z.avail_out));
            obuf = newob;
        }
        oBlockList.push(obuf);
        totalSize += obuf.length;
        if (status == Z_STREAM_END || status == Z_BUF_ERROR) {
            break;
        }
    }

    if (afterUncOffset) {
        afterUncOffset[0] = (start || 0) + z.next_in_index;
    }

    if (oBlockList.length == 1) {
        return oBlockList[0].buffer;
    } else {
        var out = new Uint8Array(totalSize);
        var cursor = 0;
        for (var i = 0; i < oBlockList.length; ++i) {
            var b = oBlockList[i];
            arrayCopy(b, 0, out, cursor, b.length);
            cursor += b.length;
        }
        return out.buffer;
    }
}

if (typeof(module) !== 'undefined') {
  module.exports = {
    inflateBuffer: jszlib_inflate_buffer,
    arrayCopy: arrayCopy
  };
}

},{}]},{},[6])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi9Vc2Vycy9kYW52ay9naXRodWIvZGFsbGlhbmNlL25vZGVfbW9kdWxlcy9ndWxwLWJyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXItcGFjay9fcHJlbHVkZS5qcyIsIi9Vc2Vycy9kYW52ay9naXRodWIvZGFsbGlhbmNlL2pzL2JhbS5qcyIsIi9Vc2Vycy9kYW52ay9naXRodWIvZGFsbGlhbmNlL2pzL2JpZ3dpZy5qcyIsIi9Vc2Vycy9kYW52ay9naXRodWIvZGFsbGlhbmNlL2pzL2Jpbi5qcyIsIi9Vc2Vycy9kYW52ay9naXRodWIvZGFsbGlhbmNlL2pzL2NvbG9yLmpzIiwiL1VzZXJzL2RhbnZrL2dpdGh1Yi9kYWxsaWFuY2UvanMvZGFzLmpzIiwiL1VzZXJzL2RhbnZrL2dpdGh1Yi9kYWxsaWFuY2UvanMvZmFrZV9lNjNmZTYyYi5qcyIsIi9Vc2Vycy9kYW52ay9naXRodWIvZGFsbGlhbmNlL2pzL2xoM3V0aWxzLmpzIiwiL1VzZXJzL2RhbnZrL2dpdGh1Yi9kYWxsaWFuY2UvanMvc2hhMS5qcyIsIi9Vc2Vycy9kYW52ay9naXRodWIvZGFsbGlhbmNlL2pzL3NwYW5zLmpzIiwiL1VzZXJzL2RhbnZrL2dpdGh1Yi9kYWxsaWFuY2UvanMvdXRpbHMuanMiLCIvVXNlcnMvZGFudmsvZ2l0aHViL2RhbGxpYW5jZS9ub2RlX21vZHVsZXMvanN6bGliL2pzL2luZmxhdGUuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pjQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqa0NBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNRQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzVIQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdjBCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3TEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUdBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ25WQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pOQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3ZEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3Rocm93IG5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIil9dmFyIGY9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGYuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sZixmLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIi8qIC0qLSBtb2RlOiBqYXZhc2NyaXB0OyBjLWJhc2ljLW9mZnNldDogNDsgaW5kZW50LXRhYnMtbW9kZTogbmlsIC0qLSAqL1xuXG4vLyBcbi8vIERhbGxpYW5jZSBHZW5vbWUgRXhwbG9yZXJcbi8vIChjKSBUaG9tYXMgRG93biAyMDA2LTIwMTFcbi8vXG4vLyBiYW0uanM6IGluZGV4ZWQgYmluYXJ5IGFsaWdubWVudHNcbi8vXG5cblwidXNlIHN0cmljdFwiO1xuXG5pZiAodHlwZW9mKHJlcXVpcmUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIHZhciBzcGFucyA9IHJlcXVpcmUoJy4vc3BhbnMnKTtcbiAgICB2YXIgUmFuZ2UgPSBzcGFucy5SYW5nZTtcbiAgICB2YXIgdW5pb24gPSBzcGFucy51bmlvbjtcbiAgICB2YXIgaW50ZXJzZWN0aW9uID0gc3BhbnMuaW50ZXJzZWN0aW9uO1xuXG4gICAgdmFyIGJpbiA9IHJlcXVpcmUoJy4vYmluJyk7XG4gICAgdmFyIHJlYWRJbnQgPSBiaW4ucmVhZEludDtcbiAgICB2YXIgcmVhZFNob3J0ID0gYmluLnJlYWRTaG9ydDtcbiAgICB2YXIgcmVhZEJ5dGUgPSBiaW4ucmVhZEJ5dGU7XG4gICAgdmFyIHJlYWRJbnQ2NCA9IGJpbi5yZWFkSW50NjQ7XG4gICAgdmFyIHJlYWRGbG9hdCA9IGJpbi5yZWFkRmxvYXQ7XG5cbiAgICB2YXIgbGgzdXRpbHMgPSByZXF1aXJlKCcuL2xoM3V0aWxzJyk7XG4gICAgdmFyIHJlYWRWb2IgPSBsaDN1dGlscy5yZWFkVm9iO1xuICAgIHZhciB1bmJnemYgPSBsaDN1dGlscy51bmJnemY7XG4gICAgdmFyIHJlZzJiaW5zID0gbGgzdXRpbHMucmVnMmJpbnM7XG4gICAgdmFyIENodW5rID0gbGgzdXRpbHMuQ2h1bms7XG59XG5cblxudmFyIEJBTV9NQUdJQyA9IDB4MTRkNDE0MjtcbnZhciBCQUlfTUFHSUMgPSAweDE0OTQxNDI7XG5cbnZhciBCYW1GbGFncyA9IHtcbiAgICBNVUxUSVBMRV9TRUdNRU5UUzogICAgICAgMHgxLFxuICAgIEFMTF9TRUdNRU5UU19BTElHTjogICAgICAweDIsXG4gICAgU0VHTUVOVF9VTk1BUFBFRDogICAgICAgIDB4NCxcbiAgICBORVhUX1NFR01FTlRfVU5NQVBQRUQ6ICAgMHg4LFxuICAgIFJFVkVSU0VfQ09NUExFTUVOVDogICAgICAweDEwLFxuICAgIE5FWFRfUkVWRVJTRV9DT01QTEVNRU5UOiAweDIwLFxuICAgIEZJUlNUX1NFR01FTlQ6ICAgICAgICAgICAweDQwLFxuICAgIExBU1RfU0VHTUVOVDogICAgICAgICAgICAweDgwLFxuICAgIFNFQ09OREFSWV9BTElHTk1FTlQ6ICAgICAweDEwMCxcbiAgICBRQ19GQUlMOiAgICAgICAgICAgICAgICAgMHgyMDAsXG4gICAgRFVQTElDQVRFOiAgICAgICAgICAgICAgIDB4NDAwLFxuICAgIFNVUFBMRU1FTlRBUlk6ICAgICAgICAgICAweDgwMFxufTtcblxuZnVuY3Rpb24gQmFtRmlsZSgpIHtcbn1cblxuZnVuY3Rpb24gbWFrZUJhbShkYXRhLCBiYWksIGNhbGxiYWNrKSB7XG4gICAgdmFyIGJhbSA9IG5ldyBCYW1GaWxlKCk7XG4gICAgYmFtLmRhdGEgPSBkYXRhO1xuICAgIGJhbS5iYWkgPSBiYWk7XG5cbiAgICBiYW0uYmFpLmZldGNoKGZ1bmN0aW9uKGhlYWRlcikgeyAgIC8vIERvIHdlIHJlYWxseSBuZWVkIHRvIGZldGNoIHRoZSB3aG9sZSB0aGluZz8gOi0oXG4gICAgICAgIGlmICghaGVhZGVyKSB7XG4gICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCwgXCJDb3VsZG4ndCBhY2Nlc3MgQkFJXCIpO1xuICAgICAgICB9XG5cbiAgICAgICAgdmFyIHVuY2JhID0gbmV3IFVpbnQ4QXJyYXkoaGVhZGVyKTtcbiAgICAgICAgdmFyIGJhaU1hZ2ljID0gcmVhZEludCh1bmNiYSwgMCk7XG4gICAgICAgIGlmIChiYWlNYWdpYyAhPSBCQUlfTUFHSUMpIHtcbiAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhudWxsLCAnTm90IGEgQkFJIGZpbGUsIG1hZ2ljPTB4JyArIGJhaU1hZ2ljLnRvU3RyaW5nKDE2KSk7XG4gICAgICAgIH1cblxuICAgICAgICB2YXIgbnJlZiA9IHJlYWRJbnQodW5jYmEsIDQpO1xuXG4gICAgICAgIGJhbS5pbmRpY2VzID0gW107XG5cbiAgICAgICAgdmFyIHAgPSA4O1xuICAgICAgICB2YXIgbWluQmxvY2tJbmRleCA9IDEwMDAwMDAwMDA7XG4gICAgICAgIGZvciAodmFyIHJlZiA9IDA7IHJlZiA8IG5yZWY7ICsrcmVmKSB7XG4gICAgICAgICAgICB2YXIgYmxvY2tTdGFydCA9IHA7XG4gICAgICAgICAgICB2YXIgbmJpbiA9IHJlYWRJbnQodW5jYmEsIHApOyBwICs9IDQ7XG4gICAgICAgICAgICBmb3IgKHZhciBiID0gMDsgYiA8IG5iaW47ICsrYikge1xuICAgICAgICAgICAgICAgIHZhciBiaW4gPSByZWFkSW50KHVuY2JhLCBwKTtcbiAgICAgICAgICAgICAgICB2YXIgbmNobmsgPSByZWFkSW50KHVuY2JhLCBwKzQpO1xuICAgICAgICAgICAgICAgIHAgKz0gOCArIChuY2huayAqIDE2KTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHZhciBuaW50diA9IHJlYWRJbnQodW5jYmEsIHApOyBwICs9IDQ7XG4gICAgICAgICAgICBcbiAgICAgICAgICAgIHZhciBxID0gcDtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgbmludHY7ICsraSkge1xuICAgICAgICAgICAgICAgIHZhciB2ID0gcmVhZFZvYih1bmNiYSwgcSk7IHEgKz0gODtcbiAgICAgICAgICAgICAgICBpZiAodikge1xuICAgICAgICAgICAgICAgICAgICB2YXIgYmkgPSB2LmJsb2NrO1xuICAgICAgICAgICAgICAgICAgICBpZiAodi5vZmZzZXQgPiAwKVxuICAgICAgICAgICAgICAgICAgICAgICAgYmkgKz0gNjU1MzY7XG5cbiAgICAgICAgICAgICAgICAgICAgaWYgKGJpIDwgbWluQmxvY2tJbmRleClcbiAgICAgICAgICAgICAgICAgICAgICAgIG1pbkJsb2NrSW5kZXggPSBiaTtcbiAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcCArPSAobmludHYgKiA4KTtcblxuXG4gICAgICAgICAgICBpZiAobmJpbiA+IDApIHtcbiAgICAgICAgICAgICAgICBiYW0uaW5kaWNlc1tyZWZdID0gbmV3IFVpbnQ4QXJyYXkoaGVhZGVyLCBibG9ja1N0YXJ0LCBwIC0gYmxvY2tTdGFydCk7XG4gICAgICAgICAgICB9ICAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgIH1cblxuICAgICAgICBiYW0uZGF0YS5zbGljZSgwLCBtaW5CbG9ja0luZGV4KS5mZXRjaChmdW5jdGlvbihyKSB7XG4gICAgICAgICAgICBpZiAoIXIpIHtcbiAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCwgXCJDb3VsZG4ndCBhY2Nlc3MgQkFNXCIpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgXG4gICAgICAgICAgICB2YXIgdW5jID0gdW5iZ3pmKHIsIHIuYnl0ZUxlbmd0aCk7XG4gICAgICAgICAgICB2YXIgdW5jYmEgPSBuZXcgVWludDhBcnJheSh1bmMpO1xuXG4gICAgICAgICAgICB2YXIgbWFnaWMgPSByZWFkSW50KHVuY2JhLCAwKTtcbiAgICAgICAgICAgIGlmIChtYWdpYyAhPSBCQU1fTUFHSUMpIHtcbiAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCwgXCJOb3QgYSBCQU0gZmlsZSwgbWFnaWM9MHhcIiArIG1hZ2ljLnRvU3RyaW5nKDE2KSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICB2YXIgaGVhZExlbiA9IHJlYWRJbnQodW5jYmEsIDQpO1xuICAgICAgICAgICAgdmFyIGhlYWRlciA9ICcnO1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBoZWFkTGVuOyArK2kpIHtcbiAgICAgICAgICAgICAgICBoZWFkZXIgKz0gU3RyaW5nLmZyb21DaGFyQ29kZSh1bmNiYVtpICsgOF0pO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICB2YXIgblJlZiA9IHJlYWRJbnQodW5jYmEsIGhlYWRMZW4gKyA4KTtcbiAgICAgICAgICAgIHZhciBwID0gaGVhZExlbiArIDEyO1xuXG4gICAgICAgICAgICBiYW0uY2hyVG9JbmRleCA9IHt9O1xuICAgICAgICAgICAgYmFtLmluZGV4VG9DaHIgPSBbXTtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgblJlZjsgKytpKSB7XG4gICAgICAgICAgICAgICAgdmFyIGxOYW1lID0gcmVhZEludCh1bmNiYSwgcCk7XG4gICAgICAgICAgICAgICAgdmFyIG5hbWUgPSAnJztcbiAgICAgICAgICAgICAgICBmb3IgKHZhciBqID0gMDsgaiA8IGxOYW1lLTE7ICsraikge1xuICAgICAgICAgICAgICAgICAgICBuYW1lICs9IFN0cmluZy5mcm9tQ2hhckNvZGUodW5jYmFbcCArIDQgKyBqXSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHZhciBsUmVmID0gcmVhZEludCh1bmNiYSwgcCArIGxOYW1lICsgNCk7XG4gICAgICAgICAgICAgICAgYmFtLmNoclRvSW5kZXhbbmFtZV0gPSBpO1xuICAgICAgICAgICAgICAgIGlmIChuYW1lLmluZGV4T2YoJ2NocicpID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgYmFtLmNoclRvSW5kZXhbbmFtZS5zdWJzdHJpbmcoMyldID0gaTtcbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICBiYW0uY2hyVG9JbmRleFsnY2hyJyArIG5hbWVdID0gaTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgYmFtLmluZGV4VG9DaHIucHVzaChuYW1lKTtcblxuICAgICAgICAgICAgICAgIHAgPSBwICsgOCArIGxOYW1lO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBpZiAoYmFtLmluZGljZXMpIHtcbiAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2soYmFtKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSk7XG4gICAgfSk7XG59XG5cblxuXG5CYW1GaWxlLnByb3RvdHlwZS5ibG9ja3NGb3JSYW5nZSA9IGZ1bmN0aW9uKHJlZklkLCBtaW4sIG1heCkge1xuICAgIHZhciBpbmRleCA9IHRoaXMuaW5kaWNlc1tyZWZJZF07XG4gICAgaWYgKCFpbmRleCkge1xuICAgICAgICByZXR1cm4gW107XG4gICAgfVxuXG4gICAgdmFyIGludEJpbnNMID0gcmVnMmJpbnMobWluLCBtYXgpO1xuICAgIHZhciBpbnRCaW5zID0gW107XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBpbnRCaW5zTC5sZW5ndGg7ICsraSkge1xuICAgICAgICBpbnRCaW5zW2ludEJpbnNMW2ldXSA9IHRydWU7XG4gICAgfVxuICAgIHZhciBsZWFmQ2h1bmtzID0gW10sIG90aGVyQ2h1bmtzID0gW107XG5cbiAgICB2YXIgbmJpbiA9IHJlYWRJbnQoaW5kZXgsIDApO1xuICAgIHZhciBwID0gNDtcbiAgICBmb3IgKHZhciBiID0gMDsgYiA8IG5iaW47ICsrYikge1xuICAgICAgICB2YXIgYmluID0gcmVhZEludChpbmRleCwgcCk7XG4gICAgICAgIHZhciBuY2huayA9IHJlYWRJbnQoaW5kZXgsIHArNCk7XG4vLyAgICAgICAgZGxvZygnYmluPScgKyBiaW4gKyAnOyBuY2huaz0nICsgbmNobmspO1xuICAgICAgICBwICs9IDg7XG4gICAgICAgIGlmIChpbnRCaW5zW2Jpbl0pIHtcbiAgICAgICAgICAgIGZvciAodmFyIGMgPSAwOyBjIDwgbmNobms7ICsrYykge1xuICAgICAgICAgICAgICAgIHZhciBjcyA9IHJlYWRWb2IoaW5kZXgsIHApO1xuICAgICAgICAgICAgICAgIHZhciBjZSA9IHJlYWRWb2IoaW5kZXgsIHAgKyA4KTtcbiAgICAgICAgICAgICAgICAoYmluIDwgNDY4MSA/IG90aGVyQ2h1bmtzIDogbGVhZkNodW5rcykucHVzaChuZXcgQ2h1bmsoY3MsIGNlKSk7XG4gICAgICAgICAgICAgICAgcCArPSAxNjtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHAgKz0gIChuY2huayAqIDE2KTtcbiAgICAgICAgfVxuICAgIH1cbi8vICAgIGRsb2coJ2xlYWZDaHVua3MgPSAnICsgbWluaUpTT05pZnkobGVhZkNodW5rcykpO1xuLy8gICAgZGxvZygnb3RoZXJDaHVua3MgPSAnICsgbWluaUpTT05pZnkob3RoZXJDaHVua3MpKTtcblxuICAgIHZhciBuaW50diA9IHJlYWRJbnQoaW5kZXgsIHApO1xuICAgIHZhciBsb3dlc3QgPSBudWxsO1xuICAgIHZhciBtaW5MaW4gPSBNYXRoLm1pbihtaW4+PjE0LCBuaW50diAtIDEpLCBtYXhMaW4gPSBNYXRoLm1pbihtYXg+PjE0LCBuaW50diAtIDEpO1xuICAgIGZvciAodmFyIGkgPSBtaW5MaW47IGkgPD0gbWF4TGluOyArK2kpIHtcbiAgICAgICAgdmFyIGxiID0gIHJlYWRWb2IoaW5kZXgsIHAgKyA0ICsgKGkgKiA4KSk7XG4gICAgICAgIGlmICghbGIpIHtcbiAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICB9XG4gICAgICAgIGlmICghbG93ZXN0IHx8IGxiLmJsb2NrIDwgbG93ZXN0LmJsb2NrIHx8IGxiLm9mZnNldCA8IGxvd2VzdC5vZmZzZXQpIHtcbiAgICAgICAgICAgIGxvd2VzdCA9IGxiO1xuICAgICAgICB9XG4gICAgfVxuICAgIC8vIGRsb2coJ0xvd2VzdCBMQiA9ICcgKyBsb3dlc3QpO1xuICAgIFxuICAgIHZhciBwcnVuZWRPdGhlckNodW5rcyA9IFtdO1xuICAgIGlmIChsb3dlc3QgIT0gbnVsbCkge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG90aGVyQ2h1bmtzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICB2YXIgY2huayA9IG90aGVyQ2h1bmtzW2ldO1xuICAgICAgICAgICAgaWYgKGNobmsubWF4di5ibG9jayA+PSBsb3dlc3QuYmxvY2sgJiYgY2huay5tYXh2Lm9mZnNldCA+PSBsb3dlc3Qub2Zmc2V0KSB7XG4gICAgICAgICAgICAgICAgcHJ1bmVkT3RoZXJDaHVua3MucHVzaChjaG5rKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cbiAgICAvLyBkbG9nKCdwcnVuZWRPdGhlckNodW5rcyA9ICcgKyBtaW5pSlNPTmlmeShwcnVuZWRPdGhlckNodW5rcykpO1xuICAgIG90aGVyQ2h1bmtzID0gcHJ1bmVkT3RoZXJDaHVua3M7XG5cbiAgICB2YXIgaW50Q2h1bmtzID0gW107XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBvdGhlckNodW5rcy5sZW5ndGg7ICsraSkge1xuICAgICAgICBpbnRDaHVua3MucHVzaChvdGhlckNodW5rc1tpXSk7XG4gICAgfVxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgbGVhZkNodW5rcy5sZW5ndGg7ICsraSkge1xuICAgICAgICBpbnRDaHVua3MucHVzaChsZWFmQ2h1bmtzW2ldKTtcbiAgICB9XG5cbiAgICBpbnRDaHVua3Muc29ydChmdW5jdGlvbihjMCwgYzEpIHtcbiAgICAgICAgdmFyIGRpZiA9IGMwLm1pbnYuYmxvY2sgLSBjMS5taW52LmJsb2NrO1xuICAgICAgICBpZiAoZGlmICE9IDApIHtcbiAgICAgICAgICAgIHJldHVybiBkaWY7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICByZXR1cm4gYzAubWludi5vZmZzZXQgLSBjMS5taW52Lm9mZnNldDtcbiAgICAgICAgfVxuICAgIH0pO1xuICAgIHZhciBtZXJnZWRDaHVua3MgPSBbXTtcbiAgICBpZiAoaW50Q2h1bmtzLmxlbmd0aCA+IDApIHtcbiAgICAgICAgdmFyIGN1ciA9IGludENodW5rc1swXTtcbiAgICAgICAgZm9yICh2YXIgaSA9IDE7IGkgPCBpbnRDaHVua3MubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgIHZhciBuYyA9IGludENodW5rc1tpXTtcbiAgICAgICAgICAgIGlmIChuYy5taW52LmJsb2NrID09IGN1ci5tYXh2LmJsb2NrIC8qICYmIG5jLm1pbnYub2Zmc2V0ID09IGN1ci5tYXh2Lm9mZnNldCAqLykgeyAvLyBubyBwb2ludCBzcGxpdHRpbmcgbWlkLWJsb2NrXG4gICAgICAgICAgICAgICAgY3VyID0gbmV3IENodW5rKGN1ci5taW52LCBuYy5tYXh2KTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgbWVyZ2VkQ2h1bmtzLnB1c2goY3VyKTtcbiAgICAgICAgICAgICAgICBjdXIgPSBuYztcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBtZXJnZWRDaHVua3MucHVzaChjdXIpO1xuICAgIH1cbiAgICAvLyBkbG9nKCdtZXJnZWRDaHVua3MgPSAnICsgbWluaUpTT05pZnkobWVyZ2VkQ2h1bmtzKSk7XG5cbiAgICByZXR1cm4gbWVyZ2VkQ2h1bmtzO1xufVxuXG5CYW1GaWxlLnByb3RvdHlwZS5mZXRjaCA9IGZ1bmN0aW9uKGNociwgbWluLCBtYXgsIGNhbGxiYWNrLCBvcHRzKSB7XG4gICAgdmFyIHRoaXNCID0gdGhpcztcbiAgICBvcHRzID0gb3B0cyB8fCB7fTtcblxuICAgIHZhciBjaHJJZCA9IHRoaXMuY2hyVG9JbmRleFtjaHJdO1xuICAgIHZhciBjaHVua3M7XG4gICAgaWYgKGNocklkID09PSB1bmRlZmluZWQpIHtcbiAgICAgICAgY2h1bmtzID0gW107XG4gICAgfSBlbHNlIHtcbiAgICAgICAgY2h1bmtzID0gdGhpcy5ibG9ja3NGb3JSYW5nZShjaHJJZCwgbWluLCBtYXgpO1xuICAgICAgICBpZiAoIWNodW5rcykge1xuICAgICAgICAgICAgY2FsbGJhY2sobnVsbCwgJ0Vycm9yIGluIGluZGV4IGZldGNoJyk7XG4gICAgICAgIH1cbiAgICB9XG4gICAgXG4gICAgdmFyIHJlY29yZHMgPSBbXTtcbiAgICB2YXIgaW5kZXggPSAwO1xuICAgIHZhciBkYXRhO1xuXG4gICAgZnVuY3Rpb24gdHJhbXAoKSB7XG4gICAgICAgIGlmIChpbmRleCA+PSBjaHVua3MubGVuZ3RoKSB7XG4gICAgICAgICAgICByZXR1cm4gY2FsbGJhY2socmVjb3Jkcyk7XG4gICAgICAgIH0gZWxzZSBpZiAoIWRhdGEpIHtcbiAgICAgICAgICAgIC8vIGRsb2coJ2ZldGNoaW5nICcgKyBpbmRleCk7XG4gICAgICAgICAgICB2YXIgYyA9IGNodW5rc1tpbmRleF07XG4gICAgICAgICAgICB2YXIgZmV0Y2hNaW4gPSBjLm1pbnYuYmxvY2s7XG4gICAgICAgICAgICB2YXIgZmV0Y2hNYXggPSBjLm1heHYuYmxvY2sgKyAoMTw8MTYpOyAvLyAqc2lnaCpcbiAgICAgICAgICAgIHRoaXNCLmRhdGEuc2xpY2UoZmV0Y2hNaW4sIGZldGNoTWF4IC0gZmV0Y2hNaW4pLmZldGNoKGZ1bmN0aW9uKHIpIHtcbiAgICAgICAgICAgICAgICBkYXRhID0gdW5iZ3pmKHIsIGMubWF4di5ibG9jayAtIGMubWludi5ibG9jayArIDEpO1xuICAgICAgICAgICAgICAgIHJldHVybiB0cmFtcCgpO1xuICAgICAgICAgICAgfSk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShkYXRhKTtcbiAgICAgICAgICAgIHRoaXNCLnJlYWRCYW1SZWNvcmRzKGJhLCBjaHVua3NbaW5kZXhdLm1pbnYub2Zmc2V0LCByZWNvcmRzLCBtaW4sIG1heCwgY2hySWQsIG9wdHMpO1xuICAgICAgICAgICAgZGF0YSA9IG51bGw7XG4gICAgICAgICAgICArK2luZGV4O1xuICAgICAgICAgICAgcmV0dXJuIHRyYW1wKCk7XG4gICAgICAgIH1cbiAgICB9XG4gICAgdHJhbXAoKTtcbn1cblxudmFyIFNFUVJFVF9ERUNPREVSID0gWyc9JywgJ0EnLCAnQycsICd4JywgJ0cnLCAneCcsICd4JywgJ3gnLCAnVCcsICd4JywgJ3gnLCAneCcsICd4JywgJ3gnLCAneCcsICdOJ107XG52YXIgQ0lHQVJfREVDT0RFUiA9IFsnTScsICdJJywgJ0QnLCAnTicsICdTJywgJ0gnLCAnUCcsICc9JywgJ1gnLCAnPycsICc/JywgJz8nLCAnPycsICc/JywgJz8nLCAnPyddO1xuXG5mdW5jdGlvbiBCYW1SZWNvcmQoKSB7XG59XG5cbkJhbUZpbGUucHJvdG90eXBlLnJlYWRCYW1SZWNvcmRzID0gZnVuY3Rpb24oYmEsIG9mZnNldCwgc2luaywgbWluLCBtYXgsIGNocklkLCBvcHRzKSB7XG4gICAgd2hpbGUgKHRydWUpIHtcbiAgICAgICAgdmFyIGJsb2NrU2l6ZSA9IHJlYWRJbnQoYmEsIG9mZnNldCk7XG4gICAgICAgIHZhciBibG9ja0VuZCA9IG9mZnNldCArIGJsb2NrU2l6ZSArIDQ7XG4gICAgICAgIGlmIChibG9ja0VuZCA+PSBiYS5sZW5ndGgpIHtcbiAgICAgICAgICAgIHJldHVybiBzaW5rO1xuICAgICAgICB9XG5cbiAgICAgICAgdmFyIHJlY29yZCA9IG5ldyBCYW1SZWNvcmQoKTtcblxuICAgICAgICB2YXIgcmVmSUQgPSByZWFkSW50KGJhLCBvZmZzZXQgKyA0KTtcbiAgICAgICAgdmFyIHBvcyA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDgpO1xuICAgICAgICBcbiAgICAgICAgdmFyIGJtbiA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDEyKTtcbiAgICAgICAgdmFyIGJpbiA9IChibW4gJiAweGZmZmYwMDAwKSA+PiAxNjtcbiAgICAgICAgdmFyIG1xID0gKGJtbiAmIDB4ZmYwMCkgPj4gODtcbiAgICAgICAgdmFyIG5sID0gYm1uICYgMHhmZjtcblxuICAgICAgICB2YXIgZmxhZ19uYyA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDE2KTtcbiAgICAgICAgdmFyIGZsYWcgPSAoZmxhZ19uYyAmIDB4ZmZmZjAwMDApID4+IDE2O1xuICAgICAgICB2YXIgbmMgPSBmbGFnX25jICYgMHhmZmZmO1xuICAgIFxuICAgICAgICB2YXIgbHNlcSA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDIwKTtcbiAgICAgICAgXG4gICAgICAgIHZhciBuZXh0UmVmICA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDI0KTtcbiAgICAgICAgdmFyIG5leHRQb3MgPSByZWFkSW50KGJhLCBvZmZzZXQgKyAyOCk7XG4gICAgICAgIFxuICAgICAgICB2YXIgdGxlbiA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDMyKTtcbiAgICBcbiAgICAgICAgcmVjb3JkLnNlZ21lbnQgPSB0aGlzLmluZGV4VG9DaHJbcmVmSURdO1xuICAgICAgICByZWNvcmQuZmxhZyA9IGZsYWc7XG4gICAgICAgIHJlY29yZC5wb3MgPSBwb3M7XG4gICAgICAgIHJlY29yZC5tcSA9IG1xO1xuICAgICAgICBpZiAob3B0cy5saWdodClcbiAgICAgICAgICAgIHJlY29yZC5zZXFMZW5ndGggPSBsc2VxO1xuXG4gICAgICAgIGlmICghb3B0cy5saWdodCkge1xuICAgICAgICAgICAgaWYgKG5leHRSZWYgPj0gMCkge1xuICAgICAgICAgICAgICAgIHJlY29yZC5uZXh0U2VnbWVudCA9IHRoaXMuaW5kZXhUb0NocltuZXh0UmVmXTtcbiAgICAgICAgICAgICAgICByZWNvcmQubmV4dFBvcyA9IG5leHRQb3M7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHZhciByZWFkTmFtZSA9ICcnO1xuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCBubC0xOyArK2opIHtcbiAgICAgICAgICAgICAgICByZWFkTmFtZSArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGJhW29mZnNldCArIDM2ICsgal0pO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmVjb3JkLnJlYWROYW1lID0gcmVhZE5hbWU7XG4gICAgICAgIFxuICAgICAgICAgICAgdmFyIHAgPSBvZmZzZXQgKyAzNiArIG5sO1xuXG4gICAgICAgICAgICB2YXIgY2lnYXIgPSAnJztcbiAgICAgICAgICAgIGZvciAodmFyIGMgPSAwOyBjIDwgbmM7ICsrYykge1xuICAgICAgICAgICAgICAgIHZhciBjaWdvcCA9IHJlYWRJbnQoYmEsIHApO1xuICAgICAgICAgICAgICAgIGNpZ2FyID0gY2lnYXIgKyAoY2lnb3A+PjQpICsgQ0lHQVJfREVDT0RFUltjaWdvcCAmIDB4Zl07XG4gICAgICAgICAgICAgICAgcCArPSA0O1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmVjb3JkLmNpZ2FyID0gY2lnYXI7XG4gICAgICAgIFxuICAgICAgICAgICAgdmFyIHNlcSA9ICcnO1xuICAgICAgICAgICAgdmFyIHNlcUJ5dGVzID0gKGxzZXEgKyAxKSA+PiAxO1xuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCBzZXFCeXRlczsgKytqKSB7XG4gICAgICAgICAgICAgICAgdmFyIHNiID0gYmFbcCArIGpdO1xuICAgICAgICAgICAgICAgIHNlcSArPSBTRVFSRVRfREVDT0RFUlsoc2IgJiAweGYwKSA+PiA0XTtcbiAgICAgICAgICAgICAgICBzZXEgKz0gU0VRUkVUX0RFQ09ERVJbKHNiICYgMHgwZildO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcCArPSBzZXFCeXRlcztcbiAgICAgICAgICAgIHJlY29yZC5zZXEgPSBzZXE7XG5cbiAgICAgICAgICAgIHZhciBxc2VxID0gJyc7XG4gICAgICAgICAgICBmb3IgKHZhciBqID0gMDsgaiA8IGxzZXE7ICsraikge1xuICAgICAgICAgICAgICAgIHFzZXEgKz0gU3RyaW5nLmZyb21DaGFyQ29kZShiYVtwICsgal0gKyAzMyk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBwICs9IGxzZXE7XG4gICAgICAgICAgICByZWNvcmQucXVhbHMgPSBxc2VxO1xuXG4gICAgICAgICAgICB3aGlsZSAocCA8IGJsb2NrRW5kKSB7XG4gICAgICAgICAgICAgICAgdmFyIHRhZyA9IFN0cmluZy5mcm9tQ2hhckNvZGUoYmFbcF0sIGJhW3AgKyAxXSk7XG4gICAgICAgICAgICAgICAgdmFyIHR5cGUgPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGJhW3AgKyAyXSk7XG4gICAgICAgICAgICAgICAgdmFyIHZhbHVlO1xuXG4gICAgICAgICAgICAgICAgaWYgKHR5cGUgPT0gJ0EnKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhbHVlID0gU3RyaW5nLmZyb21DaGFyQ29kZShiYVtwICsgM10pO1xuICAgICAgICAgICAgICAgICAgICBwICs9IDQ7XG4gICAgICAgICAgICAgICAgfSBlbHNlIGlmICh0eXBlID09ICdpJyB8fCB0eXBlID09ICdJJykge1xuICAgICAgICAgICAgICAgICAgICB2YWx1ZSA9IHJlYWRJbnQoYmEsIHAgKyAzKTtcbiAgICAgICAgICAgICAgICAgICAgcCArPSA3O1xuICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PSAnYycgfHwgdHlwZSA9PSAnQycpIHtcbiAgICAgICAgICAgICAgICAgICAgdmFsdWUgPSBiYVtwICsgM107XG4gICAgICAgICAgICAgICAgICAgIHAgKz0gNDtcbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHR5cGUgPT0gJ3MnIHx8IHR5cGUgPT0gJ1MnKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhbHVlID0gcmVhZFNob3J0KGJhLCBwICsgMyk7XG4gICAgICAgICAgICAgICAgICAgIHAgKz0gNTtcbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHR5cGUgPT0gJ2YnKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhbHVlID0gcmVhZEZsb2F0KGJhLCBwICsgMyk7XG4gICAgICAgICAgICAgICAgICAgIHAgKz0gNztcbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHR5cGUgPT0gJ1onIHx8IHR5cGUgPT0gJ0gnKSB7XG4gICAgICAgICAgICAgICAgICAgIHAgKz0gMztcbiAgICAgICAgICAgICAgICAgICAgdmFsdWUgPSAnJztcbiAgICAgICAgICAgICAgICAgICAgZm9yICg7Oykge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGNjID0gYmFbcCsrXTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChjYyA9PSAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhbHVlICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoY2MpO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfSBlbHNlIGlmICh0eXBlID09ICdCJykge1xuICAgICAgICAgICAgICAgICAgICB2YXIgYXR5cGUgPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGJhW3AgKyAzXSk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBhbGVuID0gcmVhZEludChiYSwgcCArIDQpO1xuICAgICAgICAgICAgICAgICAgICB2YXIgZWxlbjtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHJlYWRlcjtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGF0eXBlID09ICdpJyB8fCBhdHlwZSA9PSAnSScgfHwgYXR5cGUgPT0gJ2YnKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBlbGVuID0gNDtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChhdHlwZSA9PSAnZicpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgcmVhZGVyID0gcmVhZEZsb2F0O1xuICAgICAgICAgICAgICAgICAgICAgICAgZWxzZVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJlYWRlciA9IHJlYWRJbnQ7XG4gICAgICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAoYXR5cGUgPT0gJ3MnIHx8IGF0eXBlID09ICdTJykge1xuICAgICAgICAgICAgICAgICAgICAgICAgZWxlbiA9IDI7XG4gICAgICAgICAgICAgICAgICAgICAgICByZWFkZXIgPSByZWFkU2hvcnQ7XG4gICAgICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAoYXR5cGUgPT0gJ2MnIHx8IGF0eXBlID09ICdDJykge1xuICAgICAgICAgICAgICAgICAgICAgICAgZWxlbiA9IDE7XG4gICAgICAgICAgICAgICAgICAgICAgICByZWFkZXIgPSByZWFkQnl0ZTtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHRocm93ICdVbmtub3duIGFycmF5IHR5cGUgJyArIGF0eXBlO1xuICAgICAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAgICAgcCArPSA4O1xuICAgICAgICAgICAgICAgICAgICB2YWx1ZSA9IFtdO1xuICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGFsZW47ICsraSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFsdWUucHVzaChyZWFkZXIoYmEsIHApKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHAgKz0gZWxlbjtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIHRocm93ICdVbmtub3duIHR5cGUgJysgdHlwZTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgcmVjb3JkW3RhZ10gPSB2YWx1ZTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGlmICghbWluIHx8IHJlY29yZC5wb3MgPD0gbWF4ICYmIHJlY29yZC5wb3MgKyBsc2VxID49IG1pbikge1xuICAgICAgICAgICAgaWYgKGNocklkID09PSB1bmRlZmluZWQgfHwgcmVmSUQgPT0gY2hySWQpIHtcbiAgICAgICAgICAgICAgICBzaW5rLnB1c2gocmVjb3JkKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBvZmZzZXQgPSBibG9ja0VuZDtcbiAgICB9XG5cbiAgICAvLyBFeGl0cyB2aWEgdG9wIG9mIGxvb3AuXG59O1xuXG5pZiAodHlwZW9mKG1vZHVsZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgICAgIG1ha2VCYW06IG1ha2VCYW0sXG4gICAgICAgIEJBTV9NQUdJQzogQkFNX01BR0lDLFxuICAgICAgICBCQUlfTUFHSUM6IEJBSV9NQUdJQyxcbiAgICAgICAgQmFtRmxhZ3M6IEJhbUZsYWdzXG4gICAgfTtcbn0iLCIvKiAtKi0gbW9kZTogamF2YXNjcmlwdDsgYy1iYXNpYy1vZmZzZXQ6IDQ7IGluZGVudC10YWJzLW1vZGU6IG5pbCAtKi0gKi9cblxuLy8gXG4vLyBEYWxsaWFuY2UgR2Vub21lIEV4cGxvcmVyXG4vLyAoYykgVGhvbWFzIERvd24gMjAwNi0yMDEwXG4vL1xuLy8gYmlnd2lnLmpzOiBpbmRleGVkIGJpbmFyeSBXSUcgKGFuZCBCRUQpIGZpbGVzXG4vL1xuXG5cInVzZSBzdHJpY3RcIjtcblxuXG5pZiAodHlwZW9mKHJlcXVpcmUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIHZhciBzcGFucyA9IHJlcXVpcmUoJy4vc3BhbnMnKTtcbiAgICB2YXIgUmFuZ2UgPSBzcGFucy5SYW5nZTtcbiAgICB2YXIgdW5pb24gPSBzcGFucy51bmlvbjtcbiAgICB2YXIgaW50ZXJzZWN0aW9uID0gc3BhbnMuaW50ZXJzZWN0aW9uO1xuXG4gICAgdmFyIGRhcyA9IHJlcXVpcmUoJy4vZGFzJyk7XG4gICAgdmFyIERBU0ZlYXR1cmUgPSBkYXMuREFTRmVhdHVyZTtcbiAgICB2YXIgREFTR3JvdXAgPSBkYXMuREFTR3JvdXA7XG5cbiAgICB2YXIgdXRpbHMgPSByZXF1aXJlKCcuL3V0aWxzJyk7XG4gICAgdmFyIHNoYWxsb3dDb3B5ID0gdXRpbHMuc2hhbGxvd0NvcHk7XG5cbiAgICB2YXIgYmluID0gcmVxdWlyZSgnLi9iaW4nKTtcbiAgICB2YXIgcmVhZEludCA9IGJpbi5yZWFkSW50O1xuXG4gICAgdmFyIGpzemxpYiA9IHJlcXVpcmUoJ2pzemxpYicpO1xuICAgIHZhciBqc3psaWJfaW5mbGF0ZV9idWZmZXIgPSBqc3psaWIuaW5mbGF0ZUJ1ZmZlcjtcbiAgICB2YXIgYXJyYXlDb3B5ID0ganN6bGliLmFycmF5Q29weTtcbn1cblxudmFyIEJJR19XSUdfTUFHSUMgPSAweDg4OEZGQzI2O1xudmFyIEJJR19XSUdfTUFHSUNfQkUgPSAweDI2RkM4Rjg4O1xudmFyIEJJR19CRURfTUFHSUMgPSAweDg3ODlGMkVCO1xudmFyIEJJR19CRURfTUFHSUNfQkUgPSAweEVCRjI4OTg3O1xuXG5cbnZhciBCSUdfV0lHX1RZUEVfR1JBUEggPSAxO1xudmFyIEJJR19XSUdfVFlQRV9WU1RFUCA9IDI7XG52YXIgQklHX1dJR19UWVBFX0ZTVEVQID0gMztcbiAgXG52YXIgTTEgPSAyNTY7XG52YXIgTTIgPSAyNTYqMjU2O1xudmFyIE0zID0gMjU2KjI1NioyNTY7XG52YXIgTTQgPSAyNTYqMjU2KjI1NioyNTY7XG5cbnZhciBCRURfQ09MT1JfUkVHRVhQID0gbmV3IFJlZ0V4cChcIl5bMC05XSssWzAtOV0rLFswLTldK1wiKTtcblxuZnVuY3Rpb24gYndnX3JlYWRPZmZzZXQoYmEsIG8pIHtcbiAgICB2YXIgb2Zmc2V0ID0gYmFbb10gKyBiYVtvKzFdKk0xICsgYmFbbysyXSpNMiArIGJhW28rM10qTTMgKyBiYVtvKzRdKk00O1xuICAgIHJldHVybiBvZmZzZXQ7XG59XG5cbmZ1bmN0aW9uIEJpZ1dpZygpIHtcbn1cblxuQmlnV2lnLnByb3RvdHlwZS5yZWFkQ2hyb21UcmVlID0gZnVuY3Rpb24oY2FsbGJhY2spIHtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuICAgIHRoaXMuY2hyb21zVG9JRHMgPSB7fTtcbiAgICB0aGlzLmlkc1RvQ2hyb21zID0ge307XG4gICAgdGhpcy5tYXhJRCA9IDA7XG5cbiAgICB2YXIgdWRvID0gdGhpcy51bnpvb21lZERhdGFPZmZzZXQ7XG4gICAgdmFyIGViID0gKHVkbyAtIHRoaXMuY2hyb21UcmVlT2Zmc2V0KSAmIDM7XG4gICAgdWRvID0gdWRvICsgNCAtIGViO1xuXG4gICAgdGhpcy5kYXRhLnNsaWNlKHRoaXMuY2hyb21UcmVlT2Zmc2V0LCB1ZG8gLSB0aGlzLmNocm9tVHJlZU9mZnNldCkuZmV0Y2goZnVuY3Rpb24oYnB0KSB7XG4gICAgICAgIHZhciBiYSA9IG5ldyBVaW50OEFycmF5KGJwdCk7XG4gICAgICAgIHZhciBzYSA9IG5ldyBJbnQxNkFycmF5KGJwdCk7XG4gICAgICAgIHZhciBsYSA9IG5ldyBJbnQzMkFycmF5KGJwdCk7XG4gICAgICAgIHZhciBicHRNYWdpYyA9IGxhWzBdO1xuICAgICAgICB2YXIgYmxvY2tTaXplID0gbGFbMV07XG4gICAgICAgIHZhciBrZXlTaXplID0gbGFbMl07XG4gICAgICAgIHZhciB2YWxTaXplID0gbGFbM107XG4gICAgICAgIHZhciBpdGVtQ291bnQgPSBid2dfcmVhZE9mZnNldChiYSwgMTYpO1xuICAgICAgICB2YXIgcm9vdE5vZGVPZmZzZXQgPSAzMjtcblxuICAgICAgICB2YXIgYnB0UmVhZE5vZGUgPSBmdW5jdGlvbihvZmZzZXQpIHtcbiAgICAgICAgICAgIHZhciBub2RlVHlwZSA9IGJhW29mZnNldF07XG4gICAgICAgICAgICB2YXIgY250ID0gc2FbKG9mZnNldC8yKSArIDFdO1xuICAgICAgICAgICAgb2Zmc2V0ICs9IDQ7XG4gICAgICAgICAgICBmb3IgKHZhciBuID0gMDsgbiA8IGNudDsgKytuKSB7XG4gICAgICAgICAgICAgICAgaWYgKG5vZGVUeXBlID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgb2Zmc2V0ICs9IGtleVNpemU7XG4gICAgICAgICAgICAgICAgICAgIHZhciBjaGlsZE9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCBvZmZzZXQpO1xuICAgICAgICAgICAgICAgICAgICBvZmZzZXQgKz0gODtcbiAgICAgICAgICAgICAgICAgICAgY2hpbGRPZmZzZXQgLT0gdGhpc0IuY2hyb21UcmVlT2Zmc2V0O1xuICAgICAgICAgICAgICAgICAgICBicHRSZWFkTm9kZShjaGlsZE9mZnNldCk7XG4gICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGtleSA9ICcnO1xuICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBraSA9IDA7IGtpIDwga2V5U2l6ZTsgKytraSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGNoYXJDb2RlID0gYmFbb2Zmc2V0KytdO1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGNoYXJDb2RlICE9IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBrZXkgKz0gU3RyaW5nLmZyb21DaGFyQ29kZShjaGFyQ29kZSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgdmFyIGNocm9tSWQgPSAoYmFbb2Zmc2V0KzNdPDwyNCkgfCAoYmFbb2Zmc2V0KzJdPDwxNikgfCAoYmFbb2Zmc2V0KzFdPDw4KSB8IChiYVtvZmZzZXQrMF0pO1xuICAgICAgICAgICAgICAgICAgICB2YXIgY2hyb21TaXplID0gKGJhW29mZnNldCArIDddPDwyNCkgfCAoYmFbb2Zmc2V0KzZdPDwxNikgfCAoYmFbb2Zmc2V0KzVdPDw4KSB8IChiYVtvZmZzZXQrNF0pO1xuICAgICAgICAgICAgICAgICAgICBvZmZzZXQgKz0gODtcblxuICAgICAgICAgICAgICAgICAgICB0aGlzQi5jaHJvbXNUb0lEc1trZXldID0gY2hyb21JZDtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGtleS5pbmRleE9mKCdjaHInKSA9PSAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzQi5jaHJvbXNUb0lEc1trZXkuc3Vic3RyKDMpXSA9IGNocm9tSWQ7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgdGhpc0IuaWRzVG9DaHJvbXNbY2hyb21JZF0gPSBrZXk7XG4gICAgICAgICAgICAgICAgICAgIHRoaXNCLm1heElEID0gTWF0aC5tYXgodGhpc0IubWF4SUQsIGNocm9tSWQpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfTtcbiAgICAgICAgYnB0UmVhZE5vZGUocm9vdE5vZGVPZmZzZXQpO1xuXG4gICAgICAgIGNhbGxiYWNrKHRoaXNCKTtcbiAgICB9KTtcbn1cblxuZnVuY3Rpb24gQmlnV2lnVmlldyhid2csIGNpclRyZWVPZmZzZXQsIGNpclRyZWVMZW5ndGgsIGlzU3VtbWFyeSkge1xuICAgIHRoaXMuYndnID0gYndnO1xuICAgIHRoaXMuY2lyVHJlZU9mZnNldCA9IGNpclRyZWVPZmZzZXQ7XG4gICAgdGhpcy5jaXJUcmVlTGVuZ3RoID0gY2lyVHJlZUxlbmd0aDtcbiAgICB0aGlzLmlzU3VtbWFyeSA9IGlzU3VtbWFyeTtcbn1cblxuXG5cbkJpZ1dpZ1ZpZXcucHJvdG90eXBlLnJlYWRXaWdEYXRhID0gZnVuY3Rpb24oY2hyTmFtZSwgbWluLCBtYXgsIGNhbGxiYWNrKSB7XG4gICAgdmFyIGNociA9IHRoaXMuYndnLmNocm9tc1RvSURzW2Nock5hbWVdO1xuICAgIGlmIChjaHIgPT09IHVuZGVmaW5lZCkge1xuICAgICAgICAvLyBOb3QgYW4gZXJyb3IgYmVjYXVzZSBzb21lIC5id2dzIHdvbid0IGhhdmUgZGF0YSBmb3IgYWxsIGNocm9tb3NvbWVzLlxuICAgICAgICByZXR1cm4gY2FsbGJhY2soW10pO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHRoaXMucmVhZFdpZ0RhdGFCeUlkKGNociwgbWluLCBtYXgsIGNhbGxiYWNrKTtcbiAgICB9XG59XG5cbkJpZ1dpZ1ZpZXcucHJvdG90eXBlLnJlYWRXaWdEYXRhQnlJZCA9IGZ1bmN0aW9uKGNociwgbWluLCBtYXgsIGNhbGxiYWNrKSB7XG4gICAgdmFyIHRoaXNCID0gdGhpcztcbiAgICBpZiAoIXRoaXMuY2lySGVhZGVyKSB7XG4gICAgICAgIHRoaXMuYndnLmRhdGEuc2xpY2UodGhpcy5jaXJUcmVlT2Zmc2V0LCA0OCkuZmV0Y2goZnVuY3Rpb24ocmVzdWx0KSB7XG4gICAgICAgICAgICB0aGlzQi5jaXJIZWFkZXIgPSByZXN1bHQ7XG4gICAgICAgICAgICB2YXIgbGEgPSBuZXcgSW50MzJBcnJheSh0aGlzQi5jaXJIZWFkZXIpO1xuICAgICAgICAgICAgdGhpc0IuY2lyQmxvY2tTaXplID0gbGFbMV07XG4gICAgICAgICAgICB0aGlzQi5yZWFkV2lnRGF0YUJ5SWQoY2hyLCBtaW4sIG1heCwgY2FsbGJhY2spO1xuICAgICAgICB9KTtcbiAgICAgICAgcmV0dXJuO1xuICAgIH1cblxuICAgIHZhciBibG9ja3NUb0ZldGNoID0gW107XG4gICAgdmFyIG91dHN0YW5kaW5nID0gMDtcblxuICAgIHZhciBiZWZvcmVCV0cgPSBEYXRlLm5vdygpO1xuXG4gICAgdmFyIGZpbHRlciA9IGZ1bmN0aW9uKGNocm9tSWQsIGZtaW4sIGZtYXgsIHRva3MpIHtcbiAgICAgICAgcmV0dXJuICgoY2hyIDwgMCB8fCBjaHJvbUlkID09IGNocikgJiYgZm1pbiA8PSBtYXggJiYgZm1heCA+PSBtaW4pO1xuICAgIH1cblxuICAgIHZhciBjaXJGb2JSZWN1ciA9IGZ1bmN0aW9uKG9mZnNldCwgbGV2ZWwpIHtcbiAgICAgICAgaWYgKHRoaXNCLmJ3Zy5pbnN0cnVtZW50KVxuICAgICAgICAgICAgY29uc29sZS5sb2coJ2xldmVsPScgKyBsZXZlbCArICc7IG9mZnNldD0nICsgb2Zmc2V0ICsgJzsgdGltZT0nICsgKERhdGUubm93KCl8MCkpO1xuXG4gICAgICAgIG91dHN0YW5kaW5nICs9IG9mZnNldC5sZW5ndGg7XG5cbiAgICAgICAgaWYgKG9mZnNldC5sZW5ndGggPT0gMSAmJiBvZmZzZXRbMF0gLSB0aGlzQi5jaXJUcmVlT2Zmc2V0ID09IDQ4ICYmIHRoaXNCLmNhY2hlZENpclJvb3QpIHtcbiAgICAgICAgICAgIGNpckZvYlJlY3VyMih0aGlzQi5jYWNoZWRDaXJSb290LCAwLCBsZXZlbCk7XG4gICAgICAgICAgICAtLW91dHN0YW5kaW5nO1xuICAgICAgICAgICAgaWYgKG91dHN0YW5kaW5nID09IDApIHtcbiAgICAgICAgICAgICAgICB0aGlzQi5mZXRjaEZlYXR1cmVzKGZpbHRlciwgYmxvY2tzVG9GZXRjaCwgY2FsbGJhY2spO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICB9XG5cbiAgICAgICAgdmFyIG1heENpckJsb2NrU3BhbiA9IDQgKyAgKHRoaXNCLmNpckJsb2NrU2l6ZSAqIDMyKTsgICAvLyBVcHBlciBib3VuZCBvbiBzaXplLCBiYXNlZCBvbiBhIGNvbXBsZXRlbHkgZnVsbCBsZWFmIG5vZGUuXG4gICAgICAgIHZhciBzcGFucztcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBvZmZzZXQubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgIHZhciBibG9ja1NwYW4gPSBuZXcgUmFuZ2Uob2Zmc2V0W2ldLCBvZmZzZXRbaV0gKyBtYXhDaXJCbG9ja1NwYW4pO1xuICAgICAgICAgICAgc3BhbnMgPSBzcGFucyA/IHVuaW9uKHNwYW5zLCBibG9ja1NwYW4pIDogYmxvY2tTcGFuO1xuICAgICAgICB9XG4gICAgICAgIFxuICAgICAgICB2YXIgZmV0Y2hSYW5nZXMgPSBzcGFucy5yYW5nZXMoKTtcbiAgICAgICAgZm9yICh2YXIgciA9IDA7IHIgPCBmZXRjaFJhbmdlcy5sZW5ndGg7ICsrcikge1xuICAgICAgICAgICAgdmFyIGZyID0gZmV0Y2hSYW5nZXNbcl07XG4gICAgICAgICAgICBjaXJGb2JTdGFydEZldGNoKG9mZnNldCwgZnIsIGxldmVsKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIHZhciBjaXJGb2JTdGFydEZldGNoID0gZnVuY3Rpb24ob2Zmc2V0LCBmciwgbGV2ZWwsIGF0dGVtcHRzKSB7XG4gICAgICAgIHZhciBsZW5ndGggPSBmci5tYXgoKSAtIGZyLm1pbigpO1xuICAgICAgICB0aGlzQi5id2cuZGF0YS5zbGljZShmci5taW4oKSwgZnIubWF4KCkgLSBmci5taW4oKSkuZmV0Y2goZnVuY3Rpb24ocmVzdWx0QnVmZmVyKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG9mZnNldC5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgICAgIGlmIChmci5jb250YWlucyhvZmZzZXRbaV0pKSB7XG4gICAgICAgICAgICAgICAgICAgIGNpckZvYlJlY3VyMihyZXN1bHRCdWZmZXIsIG9mZnNldFtpXSAtIGZyLm1pbigpLCBsZXZlbCk7XG5cbiAgICAgICAgICAgICAgICAgICAgaWYgKG9mZnNldFtpXSAtIHRoaXNCLmNpclRyZWVPZmZzZXQgPT0gNDggJiYgb2Zmc2V0W2ldIC0gZnIubWluKCkgPT0gMClcbiAgICAgICAgICAgICAgICAgICAgICAgIHRoaXNCLmNhY2hlZENpclJvb3QgPSByZXN1bHRCdWZmZXI7XG5cbiAgICAgICAgICAgICAgICAgICAgLS1vdXRzdGFuZGluZztcbiAgICAgICAgICAgICAgICAgICAgaWYgKG91dHN0YW5kaW5nID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHRoaXNCLmZldGNoRmVhdHVyZXMoZmlsdGVyLCBibG9ja3NUb0ZldGNoLCBjYWxsYmFjayk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH0pO1xuICAgIH1cblxuICAgIHZhciBjaXJGb2JSZWN1cjIgPSBmdW5jdGlvbihjaXJCbG9ja0RhdGEsIG9mZnNldCwgbGV2ZWwpIHtcbiAgICAgICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoY2lyQmxvY2tEYXRhKTtcbiAgICAgICAgdmFyIHNhID0gbmV3IEludDE2QXJyYXkoY2lyQmxvY2tEYXRhKTtcbiAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoY2lyQmxvY2tEYXRhKTtcblxuICAgICAgICB2YXIgaXNMZWFmID0gYmFbb2Zmc2V0XTtcbiAgICAgICAgdmFyIGNudCA9IHNhW29mZnNldC8yICsgMV07XG4gICAgICAgIG9mZnNldCArPSA0O1xuXG4gICAgICAgIGlmIChpc0xlYWYgIT0gMCkge1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBjbnQ7ICsraSkge1xuICAgICAgICAgICAgICAgIHZhciBsbyA9IG9mZnNldC80O1xuICAgICAgICAgICAgICAgIHZhciBzdGFydENocm9tID0gbGFbbG9dO1xuICAgICAgICAgICAgICAgIHZhciBzdGFydEJhc2UgPSBsYVtsbyArIDFdO1xuICAgICAgICAgICAgICAgIHZhciBlbmRDaHJvbSA9IGxhW2xvICsgMl07XG4gICAgICAgICAgICAgICAgdmFyIGVuZEJhc2UgPSBsYVtsbyArIDNdO1xuICAgICAgICAgICAgICAgIHZhciBibG9ja09mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCBvZmZzZXQrMTYpO1xuICAgICAgICAgICAgICAgIHZhciBibG9ja1NpemUgPSBid2dfcmVhZE9mZnNldChiYSwgb2Zmc2V0KzI0KTtcbiAgICAgICAgICAgICAgICBpZiAoKChjaHIgPCAwIHx8IHN0YXJ0Q2hyb20gPCBjaHIpIHx8IChzdGFydENocm9tID09IGNociAmJiBzdGFydEJhc2UgPD0gbWF4KSkgJiZcbiAgICAgICAgICAgICAgICAgICAgKChjaHIgPCAwIHx8IGVuZENocm9tICAgPiBjaHIpIHx8IChlbmRDaHJvbSA9PSBjaHIgJiYgZW5kQmFzZSA+PSBtaW4pKSlcbiAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgIGJsb2Nrc1RvRmV0Y2gucHVzaCh7b2Zmc2V0OiBibG9ja09mZnNldCwgc2l6ZTogYmxvY2tTaXplfSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIG9mZnNldCArPSAzMjtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHZhciByZWN1ck9mZnNldHMgPSBbXTtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgY250OyArK2kpIHtcbiAgICAgICAgICAgICAgICB2YXIgbG8gPSBvZmZzZXQvNDtcbiAgICAgICAgICAgICAgICB2YXIgc3RhcnRDaHJvbSA9IGxhW2xvXTtcbiAgICAgICAgICAgICAgICB2YXIgc3RhcnRCYXNlID0gbGFbbG8gKyAxXTtcbiAgICAgICAgICAgICAgICB2YXIgZW5kQ2hyb20gPSBsYVtsbyArIDJdO1xuICAgICAgICAgICAgICAgIHZhciBlbmRCYXNlID0gbGFbbG8gKyAzXTtcbiAgICAgICAgICAgICAgICB2YXIgYmxvY2tPZmZzZXQgPSBid2dfcmVhZE9mZnNldChiYSwgb2Zmc2V0KzE2KTtcbiAgICAgICAgICAgICAgICBpZiAoKGNociA8IDAgfHwgc3RhcnRDaHJvbSA8IGNociB8fCAoc3RhcnRDaHJvbSA9PSBjaHIgJiYgc3RhcnRCYXNlIDw9IG1heCkpICYmXG4gICAgICAgICAgICAgICAgICAgIChjaHIgPCAwIHx8IGVuZENocm9tICAgPiBjaHIgfHwgKGVuZENocm9tID09IGNociAmJiBlbmRCYXNlID49IG1pbikpKVxuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgcmVjdXJPZmZzZXRzLnB1c2goYmxvY2tPZmZzZXQpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBvZmZzZXQgKz0gMjQ7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBpZiAocmVjdXJPZmZzZXRzLmxlbmd0aCA+IDApIHtcbiAgICAgICAgICAgICAgICBjaXJGb2JSZWN1cihyZWN1ck9mZnNldHMsIGxldmVsICsgMSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9O1xuXG4gICAgY2lyRm9iUmVjdXIoW3RoaXNCLmNpclRyZWVPZmZzZXQgKyA0OF0sIDEpO1xufVxuXG5cbkJpZ1dpZ1ZpZXcucHJvdG90eXBlLmZldGNoRmVhdHVyZXMgPSBmdW5jdGlvbihmaWx0ZXIsIGJsb2Nrc1RvRmV0Y2gsIGNhbGxiYWNrKSB7XG4gICAgdmFyIHRoaXNCID0gdGhpcztcblxuICAgIGJsb2Nrc1RvRmV0Y2guc29ydChmdW5jdGlvbihiMCwgYjEpIHtcbiAgICAgICAgcmV0dXJuIChiMC5vZmZzZXR8MCkgLSAoYjEub2Zmc2V0fDApO1xuICAgIH0pO1xuXG4gICAgaWYgKGJsb2Nrc1RvRmV0Y2gubGVuZ3RoID09IDApIHtcbiAgICAgICAgY2FsbGJhY2soW10pO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHZhciBmZWF0dXJlcyA9IFtdO1xuICAgICAgICB2YXIgY3JlYXRlRmVhdHVyZSA9IGZ1bmN0aW9uKGNociwgZm1pbiwgZm1heCwgb3B0cykge1xuICAgICAgICAgICAgaWYgKCFvcHRzKSB7XG4gICAgICAgICAgICAgICAgb3B0cyA9IHt9O1xuICAgICAgICAgICAgfVxuICAgICAgICBcbiAgICAgICAgICAgIHZhciBmID0gbmV3IERBU0ZlYXR1cmUoKTtcbiAgICAgICAgICAgIGYuX2Nocm9tSWQgPSBjaHI7XG4gICAgICAgICAgICBmLnNlZ21lbnQgPSB0aGlzQi5id2cuaWRzVG9DaHJvbXNbY2hyXTtcbiAgICAgICAgICAgIGYubWluID0gZm1pbjtcbiAgICAgICAgICAgIGYubWF4ID0gZm1heDtcbiAgICAgICAgICAgIGYudHlwZSA9ICdiaWd3aWcnO1xuICAgICAgICAgICAgXG4gICAgICAgICAgICBmb3IgKHZhciBrIGluIG9wdHMpIHtcbiAgICAgICAgICAgICAgICBmW2tdID0gb3B0c1trXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIFxuICAgICAgICAgICAgZmVhdHVyZXMucHVzaChmKTtcbiAgICAgICAgfTtcblxuICAgICAgICB2YXIgdHJhbXAgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgICAgIGlmIChibG9ja3NUb0ZldGNoLmxlbmd0aCA9PSAwKSB7XG4gICAgICAgICAgICAgICAgdmFyIGFmdGVyQldHID0gRGF0ZS5ub3coKTtcbiAgICAgICAgICAgICAgICAvLyBkbG9nKCdCV0cgZmV0Y2ggdG9vayAnICsgKGFmdGVyQldHIC0gYmVmb3JlQldHKSArICdtcycpO1xuICAgICAgICAgICAgICAgIGNhbGxiYWNrKGZlYXR1cmVzKTtcbiAgICAgICAgICAgICAgICByZXR1cm47ICAvLyBqdXN0IGluIGNhc2UuLi5cbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgdmFyIGJsb2NrID0gYmxvY2tzVG9GZXRjaFswXTtcbiAgICAgICAgICAgICAgICBpZiAoYmxvY2suZGF0YSkge1xuICAgICAgICAgICAgICAgICAgICB0aGlzQi5wYXJzZUZlYXR1cmVzKGJsb2NrLmRhdGEsIGNyZWF0ZUZlYXR1cmUsIGZpbHRlcik7XG4gICAgICAgICAgICAgICAgICAgIGJsb2Nrc1RvRmV0Y2guc3BsaWNlKDAsIDEpO1xuICAgICAgICAgICAgICAgICAgICB0cmFtcCgpO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBmZXRjaFN0YXJ0ID0gYmxvY2sub2Zmc2V0O1xuICAgICAgICAgICAgICAgICAgICB2YXIgZmV0Y2hTaXplID0gYmxvY2suc2l6ZTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGJpID0gMTtcbiAgICAgICAgICAgICAgICAgICAgd2hpbGUgKGJpIDwgYmxvY2tzVG9GZXRjaC5sZW5ndGggJiYgYmxvY2tzVG9GZXRjaFtiaV0ub2Zmc2V0ID09IChmZXRjaFN0YXJ0ICsgZmV0Y2hTaXplKSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgZmV0Y2hTaXplICs9IGJsb2Nrc1RvRmV0Y2hbYmldLnNpemU7XG4gICAgICAgICAgICAgICAgICAgICAgICArK2JpO1xuICAgICAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAgICAgdGhpc0IuYndnLmRhdGEuc2xpY2UoZmV0Y2hTdGFydCwgZmV0Y2hTaXplKS5mZXRjaChmdW5jdGlvbihyZXN1bHQpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBvZmZzZXQgPSAwO1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGJpID0gMDtcbiAgICAgICAgICAgICAgICAgICAgICAgIHdoaWxlIChvZmZzZXQgPCBmZXRjaFNpemUpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgZmIgPSBibG9ja3NUb0ZldGNoW2JpXTtcbiAgICAgICAgICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBkYXRhO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmICh0aGlzQi5id2cudW5jb21wcmVzc0J1ZlNpemUgPiAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGRhdGEgPSBqc3psaWJfaW5mbGF0ZV9idWZmZXIocmVzdWx0LCBvZmZzZXQgKyAyLCBmYi5zaXplIC0gMik7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHRtcCA9IG5ldyBVaW50OEFycmF5KGZiLnNpemUpOyAgICAvLyBGSVhNRSBpcyB0aGlzIHJlYWxseSB0aGUgYmVzdCB3ZSBjYW4gZG8/XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGFycmF5Q29weShuZXcgVWludDhBcnJheShyZXN1bHQsIG9mZnNldCwgZmIuc2l6ZSksIDAsIHRtcCwgMCwgZmIuc2l6ZSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGRhdGEgPSB0bXAuYnVmZmVyO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBmYi5kYXRhID0gZGF0YTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBvZmZzZXQgKz0gZmIuc2l6ZTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICArK2JpO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgdHJhbXAoKTtcbiAgICAgICAgICAgICAgICAgICAgfSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHRyYW1wKCk7XG4gICAgfVxufVxuXG5CaWdXaWdWaWV3LnByb3RvdHlwZS5wYXJzZUZlYXR1cmVzID0gZnVuY3Rpb24oZGF0YSwgY3JlYXRlRmVhdHVyZSwgZmlsdGVyKSB7XG4gICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoZGF0YSk7XG5cbiAgICBpZiAodGhpcy5pc1N1bW1hcnkpIHtcbiAgICAgICAgdmFyIHNhID0gbmV3IEludDE2QXJyYXkoZGF0YSk7XG4gICAgICAgIHZhciBsYSA9IG5ldyBJbnQzMkFycmF5KGRhdGEpO1xuICAgICAgICB2YXIgZmEgPSBuZXcgRmxvYXQzMkFycmF5KGRhdGEpO1xuXG4gICAgICAgIHZhciBpdGVtQ291bnQgPSBkYXRhLmJ5dGVMZW5ndGgvMzI7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgaXRlbUNvdW50OyArK2kpIHtcbiAgICAgICAgICAgIHZhciBjaHJvbUlkID0gICBsYVsoaSo4KV07XG4gICAgICAgICAgICB2YXIgc3RhcnQgPSAgICAgbGFbKGkqOCkrMV07XG4gICAgICAgICAgICB2YXIgZW5kID0gICAgICAgbGFbKGkqOCkrMl07XG4gICAgICAgICAgICB2YXIgdmFsaWRDbnQgPSAgbGFbKGkqOCkrM107XG4gICAgICAgICAgICB2YXIgbWluVmFsICAgID0gZmFbKGkqOCkrNF07XG4gICAgICAgICAgICB2YXIgbWF4VmFsICAgID0gZmFbKGkqOCkrNV07XG4gICAgICAgICAgICB2YXIgc3VtRGF0YSAgID0gZmFbKGkqOCkrNl07XG4gICAgICAgICAgICB2YXIgc3VtU3FEYXRhID0gZmFbKGkqOCkrN107XG4gICAgICAgICAgICBcbiAgICAgICAgICAgIGlmIChmaWx0ZXIoY2hyb21JZCwgc3RhcnQgKyAxLCBlbmQpKSB7XG4gICAgICAgICAgICAgICAgdmFyIHN1bW1hcnlPcHRzID0ge3R5cGU6ICdiaWd3aWcnLCBzY29yZTogc3VtRGF0YS92YWxpZENudCwgbWF4U2NvcmU6IG1heFZhbH07XG4gICAgICAgICAgICAgICAgaWYgKHRoaXMuYndnLnR5cGUgPT0gJ2JpZ2JlZCcpIHtcbiAgICAgICAgICAgICAgICAgICAgc3VtbWFyeU9wdHMudHlwZSA9ICdkZW5zaXR5JztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgY3JlYXRlRmVhdHVyZShjaHJvbUlkLCBzdGFydCArIDEsIGVuZCwgc3VtbWFyeU9wdHMpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfSBlbHNlIGlmICh0aGlzLmJ3Zy50eXBlID09ICdiaWd3aWcnKSB7XG4gICAgICAgIHZhciBzYSA9IG5ldyBJbnQxNkFycmF5KGRhdGEpO1xuICAgICAgICB2YXIgbGEgPSBuZXcgSW50MzJBcnJheShkYXRhKTtcbiAgICAgICAgdmFyIGZhID0gbmV3IEZsb2F0MzJBcnJheShkYXRhKTtcblxuICAgICAgICB2YXIgY2hyb21JZCA9IGxhWzBdO1xuICAgICAgICB2YXIgYmxvY2tTdGFydCA9IGxhWzFdO1xuICAgICAgICB2YXIgYmxvY2tFbmQgPSBsYVsyXTtcbiAgICAgICAgdmFyIGl0ZW1TdGVwID0gbGFbM107XG4gICAgICAgIHZhciBpdGVtU3BhbiA9IGxhWzRdO1xuICAgICAgICB2YXIgYmxvY2tUeXBlID0gYmFbMjBdO1xuICAgICAgICB2YXIgaXRlbUNvdW50ID0gc2FbMTFdO1xuICAgICAgICBcbiAgICAgICAgaWYgKGJsb2NrVHlwZSA9PSBCSUdfV0lHX1RZUEVfRlNURVApIHtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgaXRlbUNvdW50OyArK2kpIHtcbiAgICAgICAgICAgICAgICB2YXIgc2NvcmUgPSBmYVtpICsgNl07XG4gICAgICAgICAgICAgICAgdmFyIGZtaW4gPSBibG9ja1N0YXJ0ICsgKGkqaXRlbVN0ZXApICsgMSwgZm1heCA9IGJsb2NrU3RhcnQgKyAoaSppdGVtU3RlcCkgKyBpdGVtU3BhbjtcbiAgICAgICAgICAgICAgICBpZiAoZmlsdGVyKGNocm9tSWQsIGZtaW4sIGZtYXgpKVxuICAgICAgICAgICAgICAgICAgICBjcmVhdGVGZWF0dXJlKGNocm9tSWQsIGZtaW4sIGZtYXgsIHtzY29yZTogc2NvcmV9KTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSBlbHNlIGlmIChibG9ja1R5cGUgPT0gQklHX1dJR19UWVBFX1ZTVEVQKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGl0ZW1Db3VudDsgKytpKSB7XG4gICAgICAgICAgICAgICAgdmFyIHN0YXJ0ID0gbGFbKGkqMikgKyA2XSArIDE7XG4gICAgICAgICAgICAgICAgdmFyIGVuZCA9IHN0YXJ0ICsgaXRlbVNwYW4gLSAxO1xuICAgICAgICAgICAgICAgIHZhciBzY29yZSA9IGZhWyhpKjIpICsgN107XG4gICAgICAgICAgICAgICAgaWYgKGZpbHRlcihjaHJvbUlkLCBzdGFydCwgZW5kKSlcbiAgICAgICAgICAgICAgICAgICAgY3JlYXRlRmVhdHVyZShjaHJvbUlkLCBzdGFydCwgZW5kLCB7c2NvcmU6IHNjb3JlfSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSBpZiAoYmxvY2tUeXBlID09IEJJR19XSUdfVFlQRV9HUkFQSCkge1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBpdGVtQ291bnQ7ICsraSkge1xuICAgICAgICAgICAgICAgIHZhciBzdGFydCA9IGxhWyhpKjMpICsgNl0gKyAxO1xuICAgICAgICAgICAgICAgIHZhciBlbmQgICA9IGxhWyhpKjMpICsgN107XG4gICAgICAgICAgICAgICAgdmFyIHNjb3JlID0gZmFbKGkqMykgKyA4XTtcbiAgICAgICAgICAgICAgICBpZiAoc3RhcnQgPiBlbmQpIHtcbiAgICAgICAgICAgICAgICAgICAgc3RhcnQgPSBlbmQ7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGlmIChmaWx0ZXIoY2hyb21JZCwgc3RhcnQsIGVuZCkpXG4gICAgICAgICAgICAgICAgICAgIGNyZWF0ZUZlYXR1cmUoY2hyb21JZCwgc3RhcnQsIGVuZCwge3Njb3JlOiBzY29yZX0pO1xuICAgICAgICAgICAgfVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgY29uc29sZS5sb2coJ0N1cnJlbnRseSBub3QgaGFuZGxpbmcgYndnVHlwZT0nICsgYmxvY2tUeXBlKTtcbiAgICAgICAgfVxuICAgIH0gZWxzZSBpZiAodGhpcy5id2cudHlwZSA9PSAnYmlnYmVkJykge1xuICAgICAgICB2YXIgb2Zmc2V0ID0gMDtcbiAgICAgICAgdmFyIGRmYyA9IHRoaXMuYndnLmRlZmluZWRGaWVsZENvdW50O1xuICAgICAgICB2YXIgc2NoZW1hID0gdGhpcy5id2cuc2NoZW1hO1xuXG4gICAgICAgIHdoaWxlIChvZmZzZXQgPCBiYS5sZW5ndGgpIHtcbiAgICAgICAgICAgIHZhciBjaHJvbUlkID0gKGJhW29mZnNldCszXTw8MjQpIHwgKGJhW29mZnNldCsyXTw8MTYpIHwgKGJhW29mZnNldCsxXTw8OCkgfCAoYmFbb2Zmc2V0KzBdKTtcbiAgICAgICAgICAgIHZhciBzdGFydCA9IChiYVtvZmZzZXQrN108PDI0KSB8IChiYVtvZmZzZXQrNl08PDE2KSB8IChiYVtvZmZzZXQrNV08PDgpIHwgKGJhW29mZnNldCs0XSk7XG4gICAgICAgICAgICB2YXIgZW5kID0gKGJhW29mZnNldCsxMV08PDI0KSB8IChiYVtvZmZzZXQrMTBdPDwxNikgfCAoYmFbb2Zmc2V0KzldPDw4KSB8IChiYVtvZmZzZXQrOF0pO1xuICAgICAgICAgICAgb2Zmc2V0ICs9IDEyO1xuICAgICAgICAgICAgdmFyIHJlc3QgPSAnJztcbiAgICAgICAgICAgIHdoaWxlICh0cnVlKSB7XG4gICAgICAgICAgICAgICAgdmFyIGNoID0gYmFbb2Zmc2V0KytdO1xuICAgICAgICAgICAgICAgIGlmIChjaCAhPSAwKSB7XG4gICAgICAgICAgICAgICAgICAgIHJlc3QgKz0gU3RyaW5nLmZyb21DaGFyQ29kZShjaCk7XG4gICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICB2YXIgZmVhdHVyZU9wdHMgPSB7fTtcbiAgICAgICAgICAgIFxuICAgICAgICAgICAgdmFyIGJlZENvbHVtbnM7XG4gICAgICAgICAgICBpZiAocmVzdC5sZW5ndGggPiAwKSB7XG4gICAgICAgICAgICAgICAgYmVkQ29sdW1ucyA9IHJlc3Quc3BsaXQoJ1xcdCcpO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICBiZWRDb2x1bW5zID0gW107XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBpZiAoYmVkQ29sdW1ucy5sZW5ndGggPiAwICYmIGRmYyA+IDMpIHtcbiAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy5sYWJlbCA9IGJlZENvbHVtbnNbMF07XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBpZiAoYmVkQ29sdW1ucy5sZW5ndGggPiAxICYmIGRmYyA+IDQpIHtcbiAgICAgICAgICAgICAgICB2YXIgc2NvcmUgPSBwYXJzZUludChiZWRDb2x1bW5zWzFdKTtcbiAgICAgICAgICAgICAgICBpZiAoIWlzTmFOKHNjb3JlKSlcbiAgICAgICAgICAgICAgICAgICAgZmVhdHVyZU9wdHMuc2NvcmUgPSBzY29yZTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGlmIChiZWRDb2x1bW5zLmxlbmd0aCA+IDIgJiYgZGZjID4gNSkge1xuICAgICAgICAgICAgICAgIGZlYXR1cmVPcHRzLm9yaWVudGF0aW9uID0gYmVkQ29sdW1uc1syXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGlmIChiZWRDb2x1bW5zLmxlbmd0aCA+IDUgJiYgZGZjID4gOCkge1xuICAgICAgICAgICAgICAgIHZhciBjb2xvciA9IGJlZENvbHVtbnNbNV07XG4gICAgICAgICAgICAgICAgaWYgKEJFRF9DT0xPUl9SRUdFWFAudGVzdChjb2xvcikpIHtcbiAgICAgICAgICAgICAgICAgICAgZmVhdHVyZU9wdHMuaXRlbVJnYiA9ICdyZ2IoJyArIGNvbG9yICsgJyknO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgaWYgKGJlZENvbHVtbnMubGVuZ3RoID4gZGZjLTMgJiYgc2NoZW1hKSB7XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgY29sID0gZGZjIC0gMzsgY29sIDwgYmVkQ29sdW1ucy5sZW5ndGg7ICsrY29sKSB7XG4gICAgICAgICAgICAgICAgICAgIGZlYXR1cmVPcHRzW3NjaGVtYS5maWVsZHNbY29sKzNdLm5hbWVdID0gYmVkQ29sdW1uc1tjb2xdO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgaWYgKGZpbHRlcihjaHJvbUlkLCBzdGFydCArIDEsIGVuZCwgYmVkQ29sdW1ucykpIHtcbiAgICAgICAgICAgICAgICBpZiAoZGZjIDwgMTIpIHtcbiAgICAgICAgICAgICAgICAgICAgY3JlYXRlRmVhdHVyZShjaHJvbUlkLCBzdGFydCArIDEsIGVuZCwgZmVhdHVyZU9wdHMpO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciB0aGlja1N0YXJ0ID0gYmVkQ29sdW1uc1szXXwwO1xuICAgICAgICAgICAgICAgICAgICB2YXIgdGhpY2tFbmQgICA9IGJlZENvbHVtbnNbNF18MDtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGJsb2NrQ291bnQgPSBiZWRDb2x1bW5zWzZdfDA7XG4gICAgICAgICAgICAgICAgICAgIHZhciBibG9ja1NpemVzID0gYmVkQ29sdW1uc1s3XS5zcGxpdCgnLCcpO1xuICAgICAgICAgICAgICAgICAgICB2YXIgYmxvY2tTdGFydHMgPSBiZWRDb2x1bW5zWzhdLnNwbGl0KCcsJyk7XG5cbiAgICAgICAgICAgICAgICAgICAgaWYgKGZlYXR1cmVPcHRzLmV4b25GcmFtZXMpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBleG9uRnJhbWVzID0gZmVhdHVyZU9wdHMuZXhvbkZyYW1lcy5zcGxpdCgnLCcpO1xuICAgICAgICAgICAgICAgICAgICAgICAgZmVhdHVyZU9wdHMuZXhvbkZyYW1lcyA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgZmVhdHVyZU9wdHMudHlwZSA9ICd0cmFuc2NyaXB0J1xuICAgICAgICAgICAgICAgICAgICB2YXIgZ3JwID0gbmV3IERBU0dyb3VwKCk7XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGsgaW4gZmVhdHVyZU9wdHMpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGdycFtrXSA9IGZlYXR1cmVPcHRzW2tdO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIGdycC5pZCA9IGJlZENvbHVtbnNbMF07XG4gICAgICAgICAgICAgICAgICAgIGdycC5zZWdtZW50ID0gdGhpcy5id2cuaWRzVG9DaHJvbXNbY2hyb21JZF07XG4gICAgICAgICAgICAgICAgICAgIGdycC5taW4gPSBzdGFydCArIDE7XG4gICAgICAgICAgICAgICAgICAgIGdycC5tYXggPSBlbmQ7XG4gICAgICAgICAgICAgICAgICAgIGdycC5ub3RlcyA9IFtdO1xuICAgICAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy5ncm91cHMgPSBbZ3JwXTtcblxuICAgICAgICAgICAgICAgICAgICAvLyBNb3ZpbmcgdG93YXJkcyB1c2luZyBiaWdHZW5lUHJlZCBtb2RlbCwgYnV0IHdpbGxcbiAgICAgICAgICAgICAgICAgICAgLy8gc3RpbGwgc3VwcG9ydCBvbGQgRGFsbGlhbmNlLXN0eWxlIEJFRDEyK2dlbmUtbmFtZSBmb3IgdGhlXG4gICAgICAgICAgICAgICAgICAgIC8vIGZvcmVzZWVhYmxlIGZ1dHVyZS5cbiAgICAgICAgICAgICAgICAgICAgaWYgKGJlZENvbHVtbnMubGVuZ3RoID4gOSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGdlbmVJZCA9IGZlYXR1cmVPcHRzLmdlbmVOYW1lIHx8IGJlZENvbHVtbnNbOV07XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgZ2VuZU5hbWUgPSBnZW5lSWQ7XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoYmVkQ29sdW1ucy5sZW5ndGggPiAxMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGdlbmVOYW1lID0gYmVkQ29sdW1uc1sxMF07XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoZmVhdHVyZU9wdHMuZ2VuZU5hbWUyKVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGdlbmVOYW1lID0gZmVhdHVyZU9wdHMuZ2VuZU5hbWUyO1xuXG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgZ2cgPSBzaGFsbG93Q29weShncnApO1xuICAgICAgICAgICAgICAgICAgICAgICAgZ2cuaWQgPSBnZW5lSWQ7XG4gICAgICAgICAgICAgICAgICAgICAgICBnZy5sYWJlbCA9IGdlbmVOYW1lO1xuICAgICAgICAgICAgICAgICAgICAgICAgZ2cudHlwZSA9ICdnZW5lJztcbiAgICAgICAgICAgICAgICAgICAgICAgIGZlYXR1cmVPcHRzLmdyb3Vwcy5wdXNoKGdnKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAgICAgICAgIHZhciBzcGFuTGlzdCA9IFtdO1xuICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBiID0gMDsgYiA8IGJsb2NrQ291bnQ7ICsrYikge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGJtaW4gPSAoYmxvY2tTdGFydHNbYl18MCkgKyBzdGFydDtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBibWF4ID0gYm1pbiArIChibG9ja1NpemVzW2JdfDApO1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHNwYW4gPSBuZXcgUmFuZ2UoYm1pbiwgYm1heCk7XG4gICAgICAgICAgICAgICAgICAgICAgICBzcGFuTGlzdC5wdXNoKHNwYW4pO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIHZhciBzcGFucyA9IHVuaW9uKHNwYW5MaXN0KTtcbiAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgICAgIHZhciB0c0xpc3QgPSBzcGFucy5yYW5nZXMoKTtcbiAgICAgICAgICAgICAgICAgICAgZm9yICh2YXIgcyA9IDA7IHMgPCB0c0xpc3QubGVuZ3RoOyArK3MpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciB0cyA9IHRzTGlzdFtzXTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNyZWF0ZUZlYXR1cmUoY2hyb21JZCwgdHMubWluKCkgKyAxLCB0cy5tYXgoKSwgZmVhdHVyZU9wdHMpO1xuICAgICAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAgICAgaWYgKHRoaWNrRW5kID4gdGhpY2tTdGFydCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGNvZGluZ1JlZ2lvbiA9IChmZWF0dXJlT3B0cy5vcmllbnRhdGlvbiA9PSAnKycpID9cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBuZXcgUmFuZ2UodGhpY2tTdGFydCwgdGhpY2tFbmQgKyAzKSA6XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgbmV3IFJhbmdlKHRoaWNrU3RhcnQgLSAzLCB0aGlja0VuZCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgLy8gKy8tIDMgdG8gYWNjb3VudCBmb3Igc3RvcCBjb2RvblxuXG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgdGwgPSBpbnRlcnNlY3Rpb24oc3BhbnMsIGNvZGluZ1JlZ2lvbik7XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAodGwpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy50eXBlID0gJ3RyYW5zbGF0aW9uJztcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgdGxMaXN0ID0gdGwucmFuZ2VzKCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHJlYWRpbmdGcmFtZSA9IDA7XG5cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgdGxPZmZzZXQgPSAwO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHdoaWxlICh0bExpc3RbMF0ubWluKCkgPiB0c0xpc3RbdGxPZmZzZXRdLm1heCgpKVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0bE9mZnNldCsrO1xuXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZm9yICh2YXIgcyA9IDA7IHMgPCB0bExpc3QubGVuZ3RoOyArK3MpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgLy8gUmVjb3JkIHJlYWRpbmcgZnJhbWUgZm9yIGV2ZXJ5IGV4b25cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGluZGV4ID0gcztcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGZlYXR1cmVPcHRzLm9yaWVudGF0aW9uID09ICctJylcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluZGV4ID0gdGxMaXN0Lmxlbmd0aCAtIHMgLSAxO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgdHMgPSB0bExpc3RbaW5kZXhdO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy5yZWFkZnJhbWUgPSByZWFkaW5nRnJhbWU7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmIChleG9uRnJhbWVzKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgYnJmID0gcGFyc2VJbnQoZXhvbkZyYW1lc1tpbmRleCArIHRsT2Zmc2V0XSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAodHlwZW9mKGJyZikgPT09ICdudW1iZXInICYmIGJyZiA+PSAwICYmIGJyZiA8PSAyKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgZmVhdHVyZU9wdHMucmVhZGZyYW1lID0gYnJmO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGZlYXR1cmVPcHRzLnJlYWRmcmFtZUV4cGxpY2l0ID0gdHJ1ZTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgbGVuZ3RoID0gdHMubWF4KCkgLSB0cy5taW4oKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgcmVhZGluZ0ZyYW1lID0gKHJlYWRpbmdGcmFtZSArIGxlbmd0aCkgJSAzO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjcmVhdGVGZWF0dXJlKGNocm9tSWQsIHRzLm1pbigpICsgMSwgdHMubWF4KCksIGZlYXR1cmVPcHRzKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgICB0aHJvdyBFcnJvcihcIkRvbid0IGtub3cgd2hhdCB0byBkbyB3aXRoIFwiICsgdGhpcy5id2cudHlwZSk7XG4gICAgfVxufVxuXG4vL1xuLy8gbmFzdHkgY3V0L3Bhc3RlLCBzaG91bGQgcm9sbCBiYWNrIGluIVxuLy9cblxuQmlnV2lnVmlldy5wcm90b3R5cGUuZ2V0Rmlyc3RBZGphY2VudCA9IGZ1bmN0aW9uKGNock5hbWUsIHBvcywgZGlyLCBjYWxsYmFjaykge1xuICAgIHZhciBjaHIgPSB0aGlzLmJ3Zy5jaHJvbXNUb0lEc1tjaHJOYW1lXTtcbiAgICBpZiAoY2hyID09PSB1bmRlZmluZWQpIHtcbiAgICAgICAgLy8gTm90IGFuIGVycm9yIGJlY2F1c2Ugc29tZSAuYndncyB3b24ndCBoYXZlIGRhdGEgZm9yIGFsbCBjaHJvbW9zb21lcy5cbiAgICAgICAgcmV0dXJuIGNhbGxiYWNrKFtdKTtcbiAgICB9IGVsc2Uge1xuICAgICAgICB0aGlzLmdldEZpcnN0QWRqYWNlbnRCeUlkKGNociwgcG9zLCBkaXIsIGNhbGxiYWNrKTtcbiAgICB9XG59XG5cbkJpZ1dpZ1ZpZXcucHJvdG90eXBlLmdldEZpcnN0QWRqYWNlbnRCeUlkID0gZnVuY3Rpb24oY2hyLCBwb3MsIGRpciwgY2FsbGJhY2spIHtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuICAgIGlmICghdGhpcy5jaXJIZWFkZXIpIHtcbiAgICAgICAgdGhpcy5id2cuZGF0YS5zbGljZSh0aGlzLmNpclRyZWVPZmZzZXQsIDQ4KS5mZXRjaChmdW5jdGlvbihyZXN1bHQpIHtcbiAgICAgICAgICAgIHRoaXNCLmNpckhlYWRlciA9IHJlc3VsdDtcbiAgICAgICAgICAgIHZhciBsYSA9IG5ldyBJbnQzMkFycmF5KHRoaXNCLmNpckhlYWRlcik7XG4gICAgICAgICAgICB0aGlzQi5jaXJCbG9ja1NpemUgPSBsYVsxXTtcbiAgICAgICAgICAgIHRoaXNCLmdldEZpcnN0QWRqYWNlbnRCeUlkKGNociwgcG9zLCBkaXIsIGNhbGxiYWNrKTtcbiAgICAgICAgfSk7XG4gICAgICAgIHJldHVybjtcbiAgICB9XG5cbiAgICB2YXIgYmxvY2tUb0ZldGNoID0gbnVsbDtcbiAgICB2YXIgYmVzdEJsb2NrQ2hyID0gLTE7XG4gICAgdmFyIGJlc3RCbG9ja09mZnNldCA9IC0xO1xuXG4gICAgdmFyIG91dHN0YW5kaW5nID0gMDtcblxuICAgIHZhciBiZWZvcmVCV0cgPSBEYXRlLm5vdygpO1xuXG4gICAgdmFyIGNpckZvYlJlY3VyID0gZnVuY3Rpb24ob2Zmc2V0LCBsZXZlbCkge1xuICAgICAgICBvdXRzdGFuZGluZyArPSBvZmZzZXQubGVuZ3RoO1xuXG4gICAgICAgIHZhciBtYXhDaXJCbG9ja1NwYW4gPSA0ICsgICh0aGlzQi5jaXJCbG9ja1NpemUgKiAzMik7ICAgLy8gVXBwZXIgYm91bmQgb24gc2l6ZSwgYmFzZWQgb24gYSBjb21wbGV0ZWx5IGZ1bGwgbGVhZiBub2RlLlxuICAgICAgICB2YXIgc3BhbnM7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgb2Zmc2V0Lmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICB2YXIgYmxvY2tTcGFuID0gbmV3IFJhbmdlKG9mZnNldFtpXSwgb2Zmc2V0W2ldICsgbWF4Q2lyQmxvY2tTcGFuKTtcbiAgICAgICAgICAgIHNwYW5zID0gc3BhbnMgPyB1bmlvbihzcGFucywgYmxvY2tTcGFuKSA6IGJsb2NrU3BhbjtcbiAgICAgICAgfVxuICAgICAgICBcbiAgICAgICAgdmFyIGZldGNoUmFuZ2VzID0gc3BhbnMucmFuZ2VzKCk7XG4gICAgICAgIGZvciAodmFyIHIgPSAwOyByIDwgZmV0Y2hSYW5nZXMubGVuZ3RoOyArK3IpIHtcbiAgICAgICAgICAgIHZhciBmciA9IGZldGNoUmFuZ2VzW3JdO1xuICAgICAgICAgICAgY2lyRm9iU3RhcnRGZXRjaChvZmZzZXQsIGZyLCBsZXZlbCk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICB2YXIgY2lyRm9iU3RhcnRGZXRjaCA9IGZ1bmN0aW9uKG9mZnNldCwgZnIsIGxldmVsLCBhdHRlbXB0cykge1xuICAgICAgICB2YXIgbGVuZ3RoID0gZnIubWF4KCkgLSBmci5taW4oKTtcbiAgICAgICAgdGhpc0IuYndnLmRhdGEuc2xpY2UoZnIubWluKCksIGZyLm1heCgpIC0gZnIubWluKCkpLmZldGNoKGZ1bmN0aW9uKHJlc3VsdEJ1ZmZlcikge1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBvZmZzZXQubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgICAgICBpZiAoZnIuY29udGFpbnMob2Zmc2V0W2ldKSkge1xuICAgICAgICAgICAgICAgICAgICBjaXJGb2JSZWN1cjIocmVzdWx0QnVmZmVyLCBvZmZzZXRbaV0gLSBmci5taW4oKSwgbGV2ZWwpO1xuICAgICAgICAgICAgICAgICAgICAtLW91dHN0YW5kaW5nO1xuICAgICAgICAgICAgICAgICAgICBpZiAob3V0c3RhbmRpbmcgPT0gMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKCFibG9ja1RvRmV0Y2gpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoZGlyID4gMCAmJiAoY2hyICE9IDAgfHwgcG9zID4gMCkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIHRoaXNCLmdldEZpcnN0QWRqYWNlbnRCeUlkKDAsIDAsIGRpciwgY2FsbGJhY2spO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAoZGlyIDwgMCAmJiAoY2hyICE9IHRoaXNCLmJ3Zy5tYXhJRCB8fCBwb3MgPCAxMDAwMDAwMDAwKSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gdGhpc0IuZ2V0Rmlyc3RBZGphY2VudEJ5SWQodGhpc0IuYndnLm1heElELCAxMDAwMDAwMDAwLCBkaXIsIGNhbGxiYWNrKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKFtdKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgICAgICAgICAgICAgdGhpc0IuZmV0Y2hGZWF0dXJlcyhmdW5jdGlvbihjaHJ4LCBmbWluLCBmbWF4LCB0b2tzKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIChkaXIgPCAwICYmIChjaHJ4IDwgY2hyIHx8IGZtYXggPCBwb3MpKSB8fCAoZGlyID4gMCAmJiAoY2hyeCA+IGNociB8fCBmbWluID4gcG9zKSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9LCBbYmxvY2tUb0ZldGNoXSwgZnVuY3Rpb24oZmVhdHVyZXMpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgYmVzdEZlYXR1cmUgPSBudWxsO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBiZXN0Q2hyID0gLTE7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGJlc3RQb3MgPSAtMTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBmaSA9IDA7IGZpIDwgZmVhdHVyZXMubGVuZ3RoOyArK2ZpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBmID0gZmVhdHVyZXNbZmldO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgY2hyeCA9IGYuX2Nocm9tSWQsIGZtaW4gPSBmLm1pbiwgZm1heCA9IGYubWF4O1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoYmVzdEZlYXR1cmUgPT0gbnVsbCB8fCAoKGRpciA8IDApICYmIChjaHJ4ID4gYmVzdENociB8fCBmbWF4ID4gYmVzdFBvcykpIHx8ICgoZGlyID4gMCkgJiYgKGNocnggPCBiZXN0Q2hyIHx8IGZtaW4gPCBiZXN0UG9zKSkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJlc3RGZWF0dXJlID0gZjtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJlc3RQb3MgPSAoZGlyIDwgMCkgPyBmbWF4IDogZm1pbjtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJlc3RDaHIgPSBjaHJ4O1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGJlc3RGZWF0dXJlICE9IG51bGwpIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2soW2Jlc3RGZWF0dXJlXSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZWxzZVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2soW10pO1xuICAgICAgICAgICAgICAgICAgICAgICAgfSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH0pO1xuICAgIH1cblxuICAgIHZhciBjaXJGb2JSZWN1cjIgPSBmdW5jdGlvbihjaXJCbG9ja0RhdGEsIG9mZnNldCwgbGV2ZWwpIHtcbiAgICAgICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoY2lyQmxvY2tEYXRhKTtcbiAgICAgICAgdmFyIHNhID0gbmV3IEludDE2QXJyYXkoY2lyQmxvY2tEYXRhKTtcbiAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoY2lyQmxvY2tEYXRhKTtcblxuICAgICAgICB2YXIgaXNMZWFmID0gYmFbb2Zmc2V0XTtcbiAgICAgICAgdmFyIGNudCA9IHNhW29mZnNldC8yICsgMV07XG4gICAgICAgIG9mZnNldCArPSA0O1xuXG4gICAgICAgIGlmIChpc0xlYWYgIT0gMCkge1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBjbnQ7ICsraSkge1xuICAgICAgICAgICAgICAgIHZhciBsbyA9IG9mZnNldC80O1xuICAgICAgICAgICAgICAgIHZhciBzdGFydENocm9tID0gbGFbbG9dO1xuICAgICAgICAgICAgICAgIHZhciBzdGFydEJhc2UgPSBsYVtsbyArIDFdO1xuICAgICAgICAgICAgICAgIHZhciBlbmRDaHJvbSA9IGxhW2xvICsgMl07XG4gICAgICAgICAgICAgICAgdmFyIGVuZEJhc2UgPSBsYVtsbyArIDNdO1xuICAgICAgICAgICAgICAgIHZhciBibG9ja09mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCBvZmZzZXQrMTYpO1xuICAgICAgICAgICAgICAgIHZhciBibG9ja1NpemUgPSBid2dfcmVhZE9mZnNldChiYSwgb2Zmc2V0KzI0KTtcbiAgICAgICAgICAgICAgICBpZiAoKGRpciA8IDAgJiYgKChzdGFydENocm9tIDwgY2hyIHx8IChzdGFydENocm9tID09IGNociAmJiBzdGFydEJhc2UgPD0gcG9zKSkpKSB8fFxuICAgICAgICAgICAgICAgICAgICAoZGlyID4gMCAmJiAoKGVuZENocm9tID4gY2hyIHx8IChlbmRDaHJvbSA9PSBjaHIgJiYgZW5kQmFzZSA+PSBwb3MpKSkpKVxuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgLy8gY29uc29sZS5sb2coJ0dvdCBhbiBpbnRlcmVzdGluZyBibG9jazogc3RhcnRCYXNlPScgKyBzdGFydENocm9tICsgJzonICsgc3RhcnRCYXNlICsgJzsgZW5kQmFzZT0nICsgZW5kQ2hyb20gKyAnOicgKyBlbmRCYXNlICsgJzsgb2Zmc2V0PScgKyBibG9ja09mZnNldCArICc7IHNpemU9JyArIGJsb2NrU2l6ZSk7XG4gICAgICAgICAgICAgICAgICAgIGlmICgvX3JhbmRvbS8uZXhlYyh0aGlzQi5id2cuaWRzVG9DaHJvbXNbc3RhcnRDaHJvbV0pKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAvLyBkbG9nKCdza2lwcGluZyByYW5kb206ICcgKyB0aGlzQi5id2cuaWRzVG9DaHJvbXNbc3RhcnRDaHJvbV0pO1xuICAgICAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKGJsb2NrVG9GZXRjaCA9PSBudWxsIHx8ICgoZGlyIDwgMCkgJiYgKGVuZENocm9tID4gYmVzdEJsb2NrQ2hyIHx8IChlbmRDaHJvbSA9PSBiZXN0QmxvY2tDaHIgJiYgZW5kQmFzZSA+IGJlc3RCbG9ja09mZnNldCkpIHx8XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgKGRpciA+IDApICYmIChzdGFydENocm9tIDwgYmVzdEJsb2NrQ2hyIHx8IChzdGFydENocm9tID09IGJlc3RCbG9ja0NociAmJiBzdGFydEJhc2UgPCBiZXN0QmxvY2tPZmZzZXQpKSkpXG4gICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIC8vICAgICAgICAgICAgICAgICAgICAgICAgZGxvZygnYmVzdCBpczogc3RhcnRCYXNlPScgKyBzdGFydENocm9tICsgJzonICsgc3RhcnRCYXNlICsgJzsgZW5kQmFzZT0nICsgZW5kQ2hyb20gKyAnOicgKyBlbmRCYXNlICsgJzsgb2Zmc2V0PScgKyBibG9ja09mZnNldCArICc7IHNpemU9JyArIGJsb2NrU2l6ZSk7XG4gICAgICAgICAgICAgICAgICAgICAgICBibG9ja1RvRmV0Y2ggPSB7b2Zmc2V0OiBibG9ja09mZnNldCwgc2l6ZTogYmxvY2tTaXplfTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGJlc3RCbG9ja09mZnNldCA9IChkaXIgPCAwKSA/IGVuZEJhc2UgOiBzdGFydEJhc2U7XG4gICAgICAgICAgICAgICAgICAgICAgICBiZXN0QmxvY2tDaHIgPSAoZGlyIDwgMCkgPyBlbmRDaHJvbSA6IHN0YXJ0Q2hyb207XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgb2Zmc2V0ICs9IDMyO1xuICAgICAgICAgICAgfVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdmFyIGJlc3RSZWN1ciA9IC0xO1xuICAgICAgICAgICAgdmFyIGJlc3RQb3MgPSAtMTtcbiAgICAgICAgICAgIHZhciBiZXN0Q2hyID0gLTE7XG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGNudDsgKytpKSB7XG4gICAgICAgICAgICAgICAgdmFyIGxvID0gb2Zmc2V0LzQ7XG4gICAgICAgICAgICAgICAgdmFyIHN0YXJ0Q2hyb20gPSBsYVtsb107XG4gICAgICAgICAgICAgICAgdmFyIHN0YXJ0QmFzZSA9IGxhW2xvICsgMV07XG4gICAgICAgICAgICAgICAgdmFyIGVuZENocm9tID0gbGFbbG8gKyAyXTtcbiAgICAgICAgICAgICAgICB2YXIgZW5kQmFzZSA9IGxhW2xvICsgM107XG4gICAgICAgICAgICAgICAgdmFyIGJsb2NrT2Zmc2V0ID0gKGxhW2xvICsgNF08PDMyKSB8IChsYVtsbyArIDVdKTtcbiAgICAgICAgICAgICAgICBpZiAoKGRpciA8IDAgJiYgKChzdGFydENocm9tIDwgY2hyIHx8IChzdGFydENocm9tID09IGNociAmJiBzdGFydEJhc2UgPD0gcG9zKSkgJiZcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIChlbmRDaHJvbSAgID49IGNocikpKSB8fFxuICAgICAgICAgICAgICAgICAgICAgKGRpciA+IDAgJiYgKChlbmRDaHJvbSA+IGNociB8fCAoZW5kQ2hyb20gPT0gY2hyICYmIGVuZEJhc2UgPj0gcG9zKSkgJiZcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAoc3RhcnRDaHJvbSA8PSBjaHIpKSkpXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBpZiAoYmVzdFJlY3VyIDwgMCB8fCBlbmRCYXNlID4gYmVzdFBvcykge1xuICAgICAgICAgICAgICAgICAgICAgICAgYmVzdFJlY3VyID0gYmxvY2tPZmZzZXQ7XG4gICAgICAgICAgICAgICAgICAgICAgICBiZXN0UG9zID0gKGRpciA8IDApID8gZW5kQmFzZSA6IHN0YXJ0QmFzZTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGJlc3RDaHIgPSAoZGlyIDwgMCkgPyBlbmRDaHJvbSA6IHN0YXJ0Q2hyb207XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgb2Zmc2V0ICs9IDI0O1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgaWYgKGJlc3RSZWN1ciA+PSAwKSB7XG4gICAgICAgICAgICAgICAgY2lyRm9iUmVjdXIoW2Jlc3RSZWN1cl0sIGxldmVsICsgMSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9O1xuICAgIFxuXG4gICAgY2lyRm9iUmVjdXIoW3RoaXNCLmNpclRyZWVPZmZzZXQgKyA0OF0sIDEpO1xufVxuXG5CaWdXaWcucHJvdG90eXBlLnJlYWRXaWdEYXRhID0gZnVuY3Rpb24oY2hyTmFtZSwgbWluLCBtYXgsIGNhbGxiYWNrKSB7XG4gICAgdGhpcy5nZXRVbnpvb21lZFZpZXcoKS5yZWFkV2lnRGF0YShjaHJOYW1lLCBtaW4sIG1heCwgY2FsbGJhY2spO1xufVxuXG5CaWdXaWcucHJvdG90eXBlLmdldFVuem9vbWVkVmlldyA9IGZ1bmN0aW9uKCkge1xuICAgIGlmICghdGhpcy51bnpvb21lZFZpZXcpIHtcbiAgICAgICAgdmFyIGNpckxlbiA9IDQwMDA7XG4gICAgICAgIHZhciBuemwgPSB0aGlzLnpvb21MZXZlbHNbMF07XG4gICAgICAgIGlmIChuemwpIHtcbiAgICAgICAgICAgIGNpckxlbiA9IHRoaXMuem9vbUxldmVsc1swXS5kYXRhT2Zmc2V0IC0gdGhpcy51bnpvb21lZEluZGV4T2Zmc2V0O1xuICAgICAgICB9XG4gICAgICAgIHRoaXMudW56b29tZWRWaWV3ID0gbmV3IEJpZ1dpZ1ZpZXcodGhpcywgdGhpcy51bnpvb21lZEluZGV4T2Zmc2V0LCBjaXJMZW4sIGZhbHNlKTtcbiAgICB9XG4gICAgcmV0dXJuIHRoaXMudW56b29tZWRWaWV3O1xufVxuXG5CaWdXaWcucHJvdG90eXBlLmdldFpvb21lZFZpZXcgPSBmdW5jdGlvbih6KSB7XG4gICAgdmFyIHpoID0gdGhpcy56b29tTGV2ZWxzW3pdO1xuICAgIGlmICghemgudmlldykge1xuICAgICAgICB6aC52aWV3ID0gbmV3IEJpZ1dpZ1ZpZXcodGhpcywgemguaW5kZXhPZmZzZXQsIC8qIHRoaXMuem9vbUxldmVsc1t6ICsgMV0uZGF0YU9mZnNldCAtIHpoLmluZGV4T2Zmc2V0ICovIDQwMDAsIHRydWUpO1xuICAgIH1cbiAgICByZXR1cm4gemgudmlldztcbn1cblxuZnVuY3Rpb24gbWFrZUJ3ZyhkYXRhLCBjYWxsYmFjaywgbmFtZSkge1xuICAgIHZhciBid2cgPSBuZXcgQmlnV2lnKCk7XG4gICAgYndnLmRhdGEgPSBkYXRhO1xuICAgIGJ3Zy5uYW1lID0gbmFtZTtcbiAgICBid2cuZGF0YS5zbGljZSgwLCA1MTIpLnNhbHRlZCgpLmZldGNoKGZ1bmN0aW9uKHJlc3VsdCkge1xuICAgICAgICBpZiAoIXJlc3VsdCkge1xuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwsIFwiQ291bGRuJ3QgZmV0Y2ggZmlsZVwiKTtcbiAgICAgICAgfVxuXG4gICAgICAgIHZhciBoZWFkZXIgPSByZXN1bHQ7XG4gICAgICAgIHZhciBiYSA9IG5ldyBVaW50OEFycmF5KGhlYWRlcik7XG4gICAgICAgIHZhciBzYSA9IG5ldyBJbnQxNkFycmF5KGhlYWRlcik7XG4gICAgICAgIHZhciBsYSA9IG5ldyBJbnQzMkFycmF5KGhlYWRlcik7XG4gICAgICAgIHZhciBtYWdpYyA9IGJhWzBdICsgKE0xICogYmFbMV0pICsgKE0yICogYmFbMl0pICsgKE0zICogYmFbM10pO1xuICAgICAgICBpZiAobWFnaWMgPT0gQklHX1dJR19NQUdJQykge1xuICAgICAgICAgICAgYndnLnR5cGUgPSAnYmlnd2lnJztcbiAgICAgICAgfSBlbHNlIGlmIChtYWdpYyA9PSBCSUdfQkVEX01BR0lDKSB7XG4gICAgICAgICAgICBid2cudHlwZSA9ICdiaWdiZWQnO1xuICAgICAgICB9IGVsc2UgaWYgKG1hZ2ljID09IEJJR19XSUdfTUFHSUNfQkUgfHwgbWFnaWMgPT0gQklHX0JFRF9NQUdJQ19CRSkge1xuICAgICAgICAgICAgY2FsbGJhY2sobnVsbCwgXCJDdXJyZW50bHkgZG9uJ3Qgc3VwcG9ydCBiaWctZW5kaWFuIEJCSSBmaWxlc1wiKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGNhbGxiYWNrKG51bGwsIFwiTm90IGEgc3VwcG9ydGVkIGZvcm1hdCwgbWFnaWM9MHhcIiArIG1hZ2ljLnRvU3RyaW5nKDE2KSk7XG4gICAgICAgIH1cblxuICAgICAgICBid2cudmVyc2lvbiA9IHNhWzJdOyAgICAgICAgICAgICAvLyA0XG4gICAgICAgIGJ3Zy5udW1ab29tTGV2ZWxzID0gc2FbM107ICAgICAgIC8vIDZcbiAgICAgICAgYndnLmNocm9tVHJlZU9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCA4KTtcbiAgICAgICAgYndnLnVuem9vbWVkRGF0YU9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCAxNik7XG4gICAgICAgIGJ3Zy51bnpvb21lZEluZGV4T2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIDI0KTtcbiAgICAgICAgYndnLmZpZWxkQ291bnQgPSBzYVsxNl07ICAgICAgICAgLy8gMzJcbiAgICAgICAgYndnLmRlZmluZWRGaWVsZENvdW50ID0gc2FbMTddOyAgLy8gMzRcbiAgICAgICAgYndnLmFzT2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIDM2KTtcbiAgICAgICAgYndnLnRvdGFsU3VtbWFyeU9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCA0NCk7XG4gICAgICAgIGJ3Zy51bmNvbXByZXNzQnVmU2l6ZSA9IGxhWzEzXTsgIC8vIDUyXG4gICAgICAgIGJ3Zy5leHRIZWFkZXJPZmZzZXQgPSBid2dfcmVhZE9mZnNldChiYSwgNTYpO1xuXG4gICAgICAgIGJ3Zy56b29tTGV2ZWxzID0gW107XG4gICAgICAgIGZvciAodmFyIHpsID0gMDsgemwgPCBid2cubnVtWm9vbUxldmVsczsgKyt6bCkge1xuICAgICAgICAgICAgdmFyIHpsUmVkdWN0aW9uID0gbGFbemwqNiArIDE2XVxuICAgICAgICAgICAgdmFyIHpsRGF0YSA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCB6bCoyNCArIDcyKTtcbiAgICAgICAgICAgIHZhciB6bEluZGV4ID0gYndnX3JlYWRPZmZzZXQoYmEsIHpsKjI0ICsgODApO1xuICAgICAgICAgICAgYndnLnpvb21MZXZlbHMucHVzaCh7cmVkdWN0aW9uOiB6bFJlZHVjdGlvbiwgZGF0YU9mZnNldDogemxEYXRhLCBpbmRleE9mZnNldDogemxJbmRleH0pO1xuICAgICAgICB9XG5cbiAgICAgICAgYndnLnJlYWRDaHJvbVRyZWUoZnVuY3Rpb24oKSB7XG4gICAgICAgICAgICBid2cuZ2V0QXV0b1NRTChmdW5jdGlvbihhcykge1xuICAgICAgICAgICAgICAgIGJ3Zy5zY2hlbWEgPSBhcztcbiAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2soYndnKTtcbiAgICAgICAgICAgIH0pO1xuICAgICAgICB9KTtcbiAgICB9KTtcbn1cblxuXG5CaWdXaWcucHJvdG90eXBlLl90c0ZldGNoID0gZnVuY3Rpb24oem9vbSwgY2hyLCBtaW4sIG1heCwgY2FsbGJhY2spIHtcbiAgICB2YXIgYndnID0gdGhpcztcbiAgICBpZiAoem9vbSA+PSB0aGlzLnpvb21MZXZlbHMubGVuZ3RoIC0gMSkge1xuICAgICAgICBpZiAoIXRoaXMudG9wTGV2ZWxSZWR1Y3Rpb25DYWNoZSkge1xuICAgICAgICAgICAgdGhpcy5nZXRab29tZWRWaWV3KHRoaXMuem9vbUxldmVscy5sZW5ndGggLSAxKS5yZWFkV2lnRGF0YUJ5SWQoLTEsIDAsIDMwMDAwMDAwMCwgZnVuY3Rpb24oZmVhdHMpIHtcbiAgICAgICAgICAgICAgICBid2cudG9wTGV2ZWxSZWR1Y3Rpb25DYWNoZSA9IGZlYXRzO1xuICAgICAgICAgICAgICAgIHJldHVybiBid2cuX3RzRmV0Y2goem9vbSwgY2hyLCBtaW4sIG1heCwgY2FsbGJhY2spO1xuICAgICAgICAgICAgfSk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB2YXIgZiA9IFtdO1xuICAgICAgICAgICAgdmFyIGMgPSB0aGlzLnRvcExldmVsUmVkdWN0aW9uQ2FjaGU7XG4gICAgICAgICAgICBmb3IgKHZhciBmaSA9IDA7IGZpIDwgYy5sZW5ndGg7ICsrZmkpIHtcbiAgICAgICAgICAgICAgICBpZiAoY1tmaV0uX2Nocm9tSWQgPT0gY2hyKSB7XG4gICAgICAgICAgICAgICAgICAgIGYucHVzaChjW2ZpXSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGYpO1xuICAgICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdmFyIHZpZXc7XG4gICAgICAgIGlmICh6b29tIDwgMCkge1xuICAgICAgICAgICAgdmlldyA9IHRoaXMuZ2V0VW56b29tZWRWaWV3KCk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB2aWV3ID0gdGhpcy5nZXRab29tZWRWaWV3KHpvb20pO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiB2aWV3LnJlYWRXaWdEYXRhQnlJZChjaHIsIG1pbiwgbWF4LCBjYWxsYmFjayk7XG4gICAgfVxufVxuXG5CaWdXaWcucHJvdG90eXBlLnRocmVzaG9sZFNlYXJjaCA9IGZ1bmN0aW9uKGNock5hbWUsIHJlZmVyZW5jZVBvaW50LCBkaXIsIHRocmVzaG9sZCwgY2FsbGJhY2spIHtcbiAgICBkaXIgPSAoZGlyPDApID8gLTEgOiAxO1xuICAgIHZhciBid2cgPSB0aGlzO1xuICAgIHZhciBpbml0aWFsQ2hyID0gdGhpcy5jaHJvbXNUb0lEc1tjaHJOYW1lXTtcbiAgICB2YXIgY2FuZGlkYXRlcyA9IFt7Y2hyT3JkOiAwLCBjaHI6IGluaXRpYWxDaHIsIHpvb206IGJ3Zy56b29tTGV2ZWxzLmxlbmd0aCAtIDQsIG1pbjogMCwgbWF4OiAzMDAwMDAwMDAsIGZyb21SZWY6IHRydWV9XVxuICAgIGZvciAodmFyIGkgPSAxOyBpIDw9IHRoaXMubWF4SUQgKyAxOyArK2kpIHtcbiAgICAgICAgdmFyIGNocklkID0gKGluaXRpYWxDaHIgKyAoZGlyKmkpKSAlICh0aGlzLm1heElEICsgMSk7XG4gICAgICAgIGlmIChjaHJJZCA8IDApIFxuICAgICAgICAgICAgY2hySWQgKz0gKHRoaXMubWF4SUQgKyAxKTtcbiAgICAgICAgY2FuZGlkYXRlcy5wdXNoKHtjaHJPcmQ6IGksIGNocjogY2hySWQsIHpvb206IGJ3Zy56b29tTGV2ZWxzLmxlbmd0aCAtIDEsIG1pbjogMCwgbWF4OiAzMDAwMDAwMDB9KVxuICAgIH1cbiAgICAgICBcbiAgICBmdW5jdGlvbiBmYlRocmVzaG9sZFNlYXJjaFJlY3VyKCkge1xuICAgIFx0aWYgKGNhbmRpZGF0ZXMubGVuZ3RoID09IDApIHtcbiAgICBcdCAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgXHR9XG4gICAgXHRjYW5kaWRhdGVzLnNvcnQoZnVuY3Rpb24oYzEsIGMyKSB7XG4gICAgXHQgICAgdmFyIGQgPSBjMS56b29tIC0gYzIuem9vbTtcbiAgICBcdCAgICBpZiAoZCAhPSAwKVxuICAgIFx0XHQgICAgcmV0dXJuIGQ7XG5cbiAgICAgICAgICAgIGQgPSBjMS5jaHJPcmQgLSBjMi5jaHJPcmQ7XG4gICAgICAgICAgICBpZiAoZCAhPSAwKVxuICAgICAgICAgICAgICAgIHJldHVybiBkO1xuICAgIFx0ICAgIGVsc2VcbiAgICBcdFx0ICAgIHJldHVybiBjMS5taW4gLSBjMi5taW4gKiBkaXI7XG4gICAgXHR9KTtcblxuXHQgICAgdmFyIGNhbmRpZGF0ZSA9IGNhbmRpZGF0ZXMuc3BsaWNlKDAsIDEpWzBdO1xuICAgICAgICBid2cuX3RzRmV0Y2goY2FuZGlkYXRlLnpvb20sIGNhbmRpZGF0ZS5jaHIsIGNhbmRpZGF0ZS5taW4sIGNhbmRpZGF0ZS5tYXgsIGZ1bmN0aW9uKGZlYXRzKSB7XG4gICAgICAgICAgICB2YXIgcnAgPSBkaXIgPiAwID8gMCA6IDMwMDAwMDAwMDtcbiAgICAgICAgICAgIGlmIChjYW5kaWRhdGUuZnJvbVJlZilcbiAgICAgICAgICAgICAgICBycCA9IHJlZmVyZW5jZVBvaW50O1xuICAgICAgICAgICAgXG4gICAgICAgICAgICBmb3IgKHZhciBmaSA9IDA7IGZpIDwgZmVhdHMubGVuZ3RoOyArK2ZpKSB7XG4gICAgXHQgICAgICAgIHZhciBmID0gZmVhdHNbZmldO1xuICAgICAgICAgICAgICAgIHZhciBzY29yZTtcbiAgICAgICAgICAgICAgICBpZiAoZi5tYXhTY29yZSAhPSB1bmRlZmluZWQpXG4gICAgICAgICAgICAgICAgICAgIHNjb3JlID0gZi5tYXhTY29yZTtcbiAgICAgICAgICAgICAgICBlbHNlXG4gICAgICAgICAgICAgICAgICAgIHNjb3JlID0gZi5zY29yZTtcblxuICAgICAgICAgICAgICAgIGlmIChkaXIgPiAwKSB7XG4gICAgXHQgICAgICAgICAgICBpZiAoc2NvcmUgPiB0aHJlc2hvbGQpIHtcbiAgICAgICAgXHRcdCAgICAgICAgaWYgKGNhbmRpZGF0ZS56b29tIDwgMCkge1xuICAgICAgICBcdFx0ICAgICAgICAgICAgaWYgKGYubWluID4gcnApXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhmKTtcbiAgICAgICAgXHRcdCAgICAgICAgfSBlbHNlIGlmIChmLm1heCA+IHJwKSB7XG4gICAgICAgIFx0XHQgICAgICAgICAgICBjYW5kaWRhdGVzLnB1c2goe2NocjogY2FuZGlkYXRlLmNociwgY2hyT3JkOiBjYW5kaWRhdGUuY2hyT3JkLCB6b29tOiBjYW5kaWRhdGUuem9vbSAtIDIsIG1pbjogZi5taW4sIG1heDogZi5tYXgsIGZyb21SZWY6IGNhbmRpZGF0ZS5mcm9tUmVmfSk7XG4gICAgICAgIFx0XHQgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIGlmIChzY29yZSA+IHRocmVzaG9sZCkge1xuICAgICAgICAgICAgXHRcdCAgICBpZiAoY2FuZGlkYXRlLnpvb20gPCAwKSB7XG4gICAgICAgICAgICAgICAgXHQgICAgICAgIGlmIChmLm1heCA8IHJwKVxuICAgICAgICAgICAgICAgIFx0XHRcdCAgICByZXR1cm4gY2FsbGJhY2soZik7XG4gICAgICAgICAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKGYubWluIDwgcnApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBjYW5kaWRhdGVzLnB1c2goe2NocjogY2FuZGlkYXRlLmNociwgY2hyT3JkOiBjYW5kaWRhdGUuY2hyT3JkLCB6b29tOiBjYW5kaWRhdGUuem9vbSAtIDIsIG1pbjogZi5taW4sIG1heDogZi5tYXgsIGZyb21SZWY6IGNhbmRpZGF0ZS5mcm9tUmVmfSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgXHQgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgIFx0ICAgIH1cbiAgICAgICAgICAgIGZiVGhyZXNob2xkU2VhcmNoUmVjdXIoKTtcbiAgICAgICAgfSk7XG4gICAgfVxuICAgIFxuICAgIGZiVGhyZXNob2xkU2VhcmNoUmVjdXIoKTtcbn1cblxuQmlnV2lnLnByb3RvdHlwZS5nZXRBdXRvU1FMID0gZnVuY3Rpb24oY2FsbGJhY2spIHtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuICAgIGlmICghdGhpcy5hc09mZnNldClcbiAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwpO1xuXG5cbiAgICB0aGlzLmRhdGEuc2xpY2UodGhpcy5hc09mZnNldCwgMjA0OCkuZmV0Y2goZnVuY3Rpb24ocmVzdWx0KSB7XG4gICAgICAgIHZhciBiYSA9IG5ldyBVaW50OEFycmF5KHJlc3VsdCk7XG4gICAgICAgIHZhciBzID0gJyc7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgYmEubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgIGlmIChiYVtpXSA9PSAwKVxuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgcyArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGJhW2ldKTtcbiAgICAgICAgfVxuICAgICAgICBcbiAgICAgICAgLyogXG4gICAgICAgICAqIFF1aWNrJ24nZGlydHkgYXR0ZW1wdCB0byBwYXJzZSBhdXRvU3FsIGZvcm1hdC5cbiAgICAgICAgICogU2VlOiBodHRwOi8vd3d3LmxpbnV4am91cm5hbC5jb20vZmlsZXMvbGludXhqb3VybmFsLmNvbS9saW51eGpvdXJuYWwvYXJ0aWNsZXMvMDU5LzU5NDkvNTk0OWwyLmh0bWxcbiAgICAgICAgICovXG5cbiAgICAgICAgdmFyIGhlYWRlcl9yZSA9IC8oXFx3KylcXHMrKFxcdyspXFxzKyhcIihbXlwiXSspXCIpP1xccytcXChcXHMqLztcbiAgICAgICAgdmFyIGZpZWxkX3JlID0gLyhbXFx3XFxbXFxdXSspXFxzKyhcXHcrKVxccyo7XFxzKihcIihbXlwiXSspXCIpP1xccyovZztcblxuICAgICAgICB2YXIgaGVhZGVyTWF0Y2ggPSBoZWFkZXJfcmUuZXhlYyhzKTtcbiAgICAgICAgaWYgKGhlYWRlck1hdGNoKSB7XG4gICAgICAgICAgICB2YXIgYXMgPSB7XG4gICAgICAgICAgICAgICAgZGVjbFR5cGU6IGhlYWRlck1hdGNoWzFdLFxuICAgICAgICAgICAgICAgIG5hbWU6IGhlYWRlck1hdGNoWzJdLFxuICAgICAgICAgICAgICAgIGNvbW1lbnQ6IGhlYWRlck1hdGNoWzRdLFxuXG4gICAgICAgICAgICAgICAgZmllbGRzOiBbXVxuICAgICAgICAgICAgfTtcblxuICAgICAgICAgICAgcyA9IHMuc3Vic3RyaW5nKGhlYWRlck1hdGNoWzBdKTtcbiAgICAgICAgICAgIGZvciAodmFyIG0gPSBmaWVsZF9yZS5leGVjKHMpOyBtICE9IG51bGw7IG0gPSBmaWVsZF9yZS5leGVjKHMpKSB7XG4gICAgICAgICAgICAgICAgYXMuZmllbGRzLnB1c2goe3R5cGU6IG1bMV0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgIG5hbWU6IG1bMl0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgIGNvbW1lbnQ6IG1bNF19KTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGFzKTtcbiAgICAgICAgfVxuICAgIH0pO1xufVxuXG5CaWdXaWcucHJvdG90eXBlLmdldEV4dHJhSW5kaWNlcyA9IGZ1bmN0aW9uKGNhbGxiYWNrKSB7XG4gICAgdmFyIHRoaXNCID0gdGhpcztcbiAgICBpZiAodGhpcy52ZXJzaW9uIDwgNCB8fCB0aGlzLmV4dEhlYWRlck9mZnNldCA9PSAwIHx8IHRoaXMudHlwZSAhPSAnYmlnYmVkJykge1xuICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdGhpcy5kYXRhLnNsaWNlKHRoaXMuZXh0SGVhZGVyT2Zmc2V0LCA2NCkuZmV0Y2goZnVuY3Rpb24ocmVzdWx0KSB7XG4gICAgICAgICAgICBpZiAoIXJlc3VsdCkge1xuICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhudWxsLCBcIkNvdWxkbid0IGZldGNoIGV4dGVuc2lvbiBoZWFkZXJcIik7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHZhciBiYSA9IG5ldyBVaW50OEFycmF5KHJlc3VsdCk7XG4gICAgICAgICAgICB2YXIgc2EgPSBuZXcgSW50MTZBcnJheShyZXN1bHQpO1xuICAgICAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkocmVzdWx0KTtcbiAgICAgICAgICAgIFxuICAgICAgICAgICAgdmFyIGV4dEhlYWRlclNpemUgPSBzYVswXTtcbiAgICAgICAgICAgIHZhciBleHRyYUluZGV4Q291bnQgPSBzYVsxXTtcbiAgICAgICAgICAgIHZhciBleHRyYUluZGV4TGlzdE9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCA0KTtcblxuICAgICAgICAgICAgaWYgKGV4dHJhSW5kZXhDb3VudCA9PSAwKSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwpO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAvLyBGSVhNRSAyMGJ5dGUgcmVjb3JkcyBvbmx5IG1ha2Ugc2Vuc2UgZm9yIHNpbmdsZS1maWVsZCBpbmRpY2VzLlxuICAgICAgICAgICAgLy8gUmlnaHQgbm93LCB0aGVzZSBzZWVtIHRvIGJlIHRoZSBvbmx5IHRoaW5ncyBhcm91bmQsIGJ1dCB0aGUgZm9ybWF0XG4gICAgICAgICAgICAvLyBpcyBhY3R1YWxseSBtb3JlIGdlbmVyYWwuXG4gICAgICAgICAgICB0aGlzQi5kYXRhLnNsaWNlKGV4dHJhSW5kZXhMaXN0T2Zmc2V0LCBleHRyYUluZGV4Q291bnQgKiAyMCkuZmV0Y2goZnVuY3Rpb24oZWlsKSB7XG4gICAgICAgICAgICAgICAgaWYgKCFlaWwpIHtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwsIFwiQ291bGRuJ3QgZmV0Y2ggaW5kZXggaW5mb1wiKTtcbiAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShlaWwpO1xuICAgICAgICAgICAgICAgIHZhciBzYSA9IG5ldyBJbnQxNkFycmF5KGVpbCk7XG4gICAgICAgICAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoZWlsKTtcblxuICAgICAgICAgICAgICAgIHZhciBpbmRpY2VzID0gW107XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgaWkgPSAwOyBpaSA8IGV4dHJhSW5kZXhDb3VudDsgKytpaSkge1xuICAgICAgICAgICAgICAgICAgICB2YXIgZWlUeXBlID0gc2FbaWkqMTBdO1xuICAgICAgICAgICAgICAgICAgICB2YXIgZWlGaWVsZENvdW50ID0gc2FbaWkqMTAgKyAxXTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGVpT2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIGlpKjIwICsgNCk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBlaUZpZWxkID0gc2FbaWkqMTAgKyA4XVxuICAgICAgICAgICAgICAgICAgICB2YXIgaW5kZXggPSBuZXcgQkJJRXh0cmFJbmRleCh0aGlzQiwgZWlUeXBlLCBlaUZpZWxkQ291bnQsIGVpT2Zmc2V0LCBlaUZpZWxkKTtcbiAgICAgICAgICAgICAgICAgICAgaW5kaWNlcy5wdXNoKGluZGV4KTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgY2FsbGJhY2soaW5kaWNlcyk7XG4gICAgICAgICAgICB9KTtcbiAgICAgICAgfSk7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBCQklFeHRyYUluZGV4KGJiaSwgdHlwZSwgZmllbGRDb3VudCwgb2Zmc2V0LCBmaWVsZCkge1xuICAgIHRoaXMuYmJpID0gYmJpO1xuICAgIHRoaXMudHlwZSA9IHR5cGU7XG4gICAgdGhpcy5maWVsZENvdW50ID0gZmllbGRDb3VudDtcbiAgICB0aGlzLm9mZnNldCA9IG9mZnNldDtcbiAgICB0aGlzLmZpZWxkID0gZmllbGQ7XG59XG5cbkJCSUV4dHJhSW5kZXgucHJvdG90eXBlLmxvb2t1cCA9IGZ1bmN0aW9uKG5hbWUsIGNhbGxiYWNrKSB7XG4gICAgdmFyIHRoaXNCID0gdGhpcztcblxuICAgIHRoaXMuYmJpLmRhdGEuc2xpY2UodGhpcy5vZmZzZXQsIDMyKS5mZXRjaChmdW5jdGlvbihicHQpIHtcbiAgICAgICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIHNhID0gbmV3IEludDE2QXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIGJwdE1hZ2ljID0gbGFbMF07XG4gICAgICAgIHZhciBibG9ja1NpemUgPSBsYVsxXTtcbiAgICAgICAgdmFyIGtleVNpemUgPSBsYVsyXTtcbiAgICAgICAgdmFyIHZhbFNpemUgPSBsYVszXTtcbiAgICAgICAgdmFyIGl0ZW1Db3VudCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCAxNik7XG4gICAgICAgIHZhciByb290Tm9kZU9mZnNldCA9IDMyO1xuXG4gICAgICAgIGZ1bmN0aW9uIGJwdFJlYWROb2RlKG5vZGVPZmZzZXQpIHtcbiAgICAgICAgICAgIHRoaXNCLmJiaS5kYXRhLnNsaWNlKG5vZGVPZmZzZXQsIDQgKyAoYmxvY2tTaXplICogKGtleVNpemUgKyB2YWxTaXplKSkpLmZldGNoKGZ1bmN0aW9uKG5vZGUpIHtcbiAgICAgICAgICAgICAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShub2RlKTtcbiAgICAgICAgICAgICAgICB2YXIgc2EgPSBuZXcgVWludDE2QXJyYXkobm9kZSk7XG4gICAgICAgICAgICAgICAgdmFyIGxhID0gbmV3IFVpbnQzMkFycmF5KG5vZGUpO1xuXG4gICAgICAgICAgICAgICAgdmFyIG5vZGVUeXBlID0gYmFbMF07XG4gICAgICAgICAgICAgICAgdmFyIGNudCA9IHNhWzFdO1xuXG4gICAgICAgICAgICAgICAgdmFyIG9mZnNldCA9IDQ7XG4gICAgICAgICAgICAgICAgaWYgKG5vZGVUeXBlID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGxhc3RDaGlsZE9mZnNldCA9IG51bGw7XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIG4gPSAwOyBuIDwgY250OyArK24pIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBrZXkgPSAnJztcbiAgICAgICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGtpID0gMDsga2kgPCBrZXlTaXplOyArK2tpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGNoYXJDb2RlID0gYmFbb2Zmc2V0KytdO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmIChjaGFyQ29kZSAhPSAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGtleSArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGNoYXJDb2RlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBjaGlsZE9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCBvZmZzZXQpO1xuICAgICAgICAgICAgICAgICAgICAgICAgb2Zmc2V0ICs9IDg7XG4gICAgICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChuYW1lLmxvY2FsZUNvbXBhcmUoa2V5KSA8IDAgJiYgbGFzdENoaWxkT2Zmc2V0KSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnB0UmVhZE5vZGUobGFzdENoaWxkT2Zmc2V0KTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm47XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICBsYXN0Q2hpbGRPZmZzZXQgPSBjaGlsZE9mZnNldDtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICBicHRSZWFkTm9kZShsYXN0Q2hpbGRPZmZzZXQpO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIG4gPSAwOyBuIDwgY250OyArK24pIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBrZXkgPSAnJztcbiAgICAgICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGtpID0gMDsga2kgPCBrZXlTaXplOyArK2tpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGNoYXJDb2RlID0gYmFbb2Zmc2V0KytdO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmIChjaGFyQ29kZSAhPSAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGtleSArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGNoYXJDb2RlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgICAgIC8vIFNwZWNpZmljIGZvciBFSSBjYXNlLlxuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGtleSA9PSBuYW1lKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHN0YXJ0ID0gYndnX3JlYWRPZmZzZXQoYmEsIG9mZnNldCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGxlbmd0aCA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDgpO1xuXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIHRoaXNCLmJiaS5nZXRVbnpvb21lZFZpZXcoKS5mZXRjaEZlYXR1cmVzKFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBmdW5jdGlvbihjaHIsIG1pbiwgbWF4LCB0b2tzKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAodG9rcyAmJiB0b2tzLmxlbmd0aCA+IHRoaXNCLmZpZWxkIC0gMylcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gdG9rc1t0aGlzQi5maWVsZCAtIDNdID09IG5hbWU7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBbe29mZnNldDogc3RhcnQsIHNpemU6IGxlbmd0aH1dLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY2FsbGJhY2spO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgb2Zmc2V0ICs9IHZhbFNpemU7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKFtdKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9KTtcbiAgICAgICAgfVxuXG4gICAgICAgIGJwdFJlYWROb2RlKHRoaXNCLm9mZnNldCArIHJvb3ROb2RlT2Zmc2V0KTtcbiAgICB9KTtcbn1cblxuaWYgKHR5cGVvZihtb2R1bGUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIG1vZHVsZS5leHBvcnRzID0ge1xuICAgICAgICBtYWtlQndnOiBtYWtlQndnLFxuICAgICAgICBCSUdfQkVEX01BR0lDOiBCSUdfQkVEX01BR0lDLFxuICAgICAgICBCSUdfV0lHX01BR0lDOiBCSUdfV0lHX01BR0lDXG4gICAgfVxufVxuIiwiLyogLSotIG1vZGU6IGphdmFzY3JpcHQ7IGMtYmFzaWMtb2Zmc2V0OiA0OyBpbmRlbnQtdGFicy1tb2RlOiBuaWwgLSotICovXG5cbi8vIFxuLy8gRGFsbGlhbmNlIEdlbm9tZSBFeHBsb3JlclxuLy8gKGMpIFRob21hcyBEb3duIDIwMDYtMjAxMVxuLy9cbi8vIGJpbi5qcyBnZW5lcmFsIGJpbmFyeSBkYXRhIHN1cHBvcnRcbi8vXG5cblwidXNlIHN0cmljdFwiO1xuXG5pZiAodHlwZW9mKHJlcXVpcmUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIHZhciB1dGlscyA9IHJlcXVpcmUoJy4vdXRpbHMnKTtcbiAgICB2YXIgc2hhbGxvd0NvcHkgPSB1dGlscy5zaGFsbG93Q29weTtcblxuICAgIHZhciBzaGExID0gcmVxdWlyZSgnLi9zaGExJyk7XG4gICAgdmFyIGI2NF9zaGExID0gc2hhMS5iNjRfc2hhMTtcbn1cblxuZnVuY3Rpb24gQmxvYkZldGNoYWJsZShiKSB7XG4gICAgdGhpcy5ibG9iID0gYjtcbn1cblxuQmxvYkZldGNoYWJsZS5wcm90b3R5cGUuc2xpY2UgPSBmdW5jdGlvbihzdGFydCwgbGVuZ3RoKSB7XG4gICAgdmFyIGI7XG5cbiAgICBpZiAodGhpcy5ibG9iLnNsaWNlKSB7XG4gICAgICAgIGlmIChsZW5ndGgpIHtcbiAgICAgICAgICAgIGIgPSB0aGlzLmJsb2Iuc2xpY2Uoc3RhcnQsIHN0YXJ0ICsgbGVuZ3RoKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGIgPSB0aGlzLmJsb2Iuc2xpY2Uoc3RhcnQpO1xuICAgICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgICAgaWYgKGxlbmd0aCkge1xuICAgICAgICAgICAgYiA9IHRoaXMuYmxvYi53ZWJraXRTbGljZShzdGFydCwgc3RhcnQgKyBsZW5ndGgpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgYiA9IHRoaXMuYmxvYi53ZWJraXRTbGljZShzdGFydCk7XG4gICAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuIG5ldyBCbG9iRmV0Y2hhYmxlKGIpO1xufVxuXG5CbG9iRmV0Y2hhYmxlLnByb3RvdHlwZS5zYWx0ZWQgPSBmdW5jdGlvbigpIHtyZXR1cm4gdGhpczt9XG5cbmlmICh0eXBlb2YoRmlsZVJlYWRlcikgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgLy8gY29uc29sZS5sb2coJ2RlZmluaW5nIGFzeW5jIEJsb2JGZXRjaGFibGUuZmV0Y2gnKTtcblxuICAgIEJsb2JGZXRjaGFibGUucHJvdG90eXBlLmZldGNoID0gZnVuY3Rpb24oY2FsbGJhY2spIHtcbiAgICAgICAgdmFyIHJlYWRlciA9IG5ldyBGaWxlUmVhZGVyKCk7XG4gICAgICAgIHJlYWRlci5vbmxvYWRlbmQgPSBmdW5jdGlvbihldikge1xuICAgICAgICAgICAgY2FsbGJhY2soYnN0cmluZ1RvQnVmZmVyKHJlYWRlci5yZXN1bHQpKTtcbiAgICAgICAgfTtcbiAgICAgICAgcmVhZGVyLnJlYWRBc0JpbmFyeVN0cmluZyh0aGlzLmJsb2IpO1xuICAgIH1cblxufSBlbHNlIHtcbiAgICAvLyBpZiAoY29uc29sZSAmJiBjb25zb2xlLmxvZylcbiAgICAvLyAgICBjb25zb2xlLmxvZygnZGVmaW5pbmcgc3luYyBCbG9iRmV0Y2hhYmxlLmZldGNoJyk7XG5cbiAgICBCbG9iRmV0Y2hhYmxlLnByb3RvdHlwZS5mZXRjaCA9IGZ1bmN0aW9uKGNhbGxiYWNrKSB7XG4gICAgICAgIHZhciByZWFkZXIgPSBuZXcgRmlsZVJlYWRlclN5bmMoKTtcbiAgICAgICAgdHJ5IHtcbiAgICAgICAgICAgIHZhciByZXMgPSByZWFkZXIucmVhZEFzQXJyYXlCdWZmZXIodGhpcy5ibG9iKTtcbiAgICAgICAgICAgIGNhbGxiYWNrKHJlcyk7XG4gICAgICAgIH0gY2F0Y2ggKGUpIHtcbiAgICAgICAgICAgIGNhbGxiYWNrKG51bGwsIGUpO1xuICAgICAgICB9XG4gICAgfVxufVxuXG5mdW5jdGlvbiBVUkxGZXRjaGFibGUodXJsLCBzdGFydCwgZW5kLCBvcHRzKSB7XG4gICAgaWYgKCFvcHRzKSB7XG4gICAgICAgIGlmICh0eXBlb2Ygc3RhcnQgPT09ICdvYmplY3QnKSB7XG4gICAgICAgICAgICBvcHRzID0gc3RhcnQ7XG4gICAgICAgICAgICBzdGFydCA9IHVuZGVmaW5lZDtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIG9wdHMgPSB7fTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIHRoaXMudXJsID0gdXJsO1xuICAgIHRoaXMuc3RhcnQgPSBzdGFydCB8fCAwO1xuICAgIGlmIChlbmQpIHtcbiAgICAgICAgdGhpcy5lbmQgPSBlbmQ7XG4gICAgfVxuICAgIHRoaXMub3B0cyA9IG9wdHM7XG59XG5cblVSTEZldGNoYWJsZS5wcm90b3R5cGUuc2xpY2UgPSBmdW5jdGlvbihzLCBsKSB7XG4gICAgaWYgKHMgPCAwKSB7XG4gICAgICAgIHRocm93ICdCYWQgc2xpY2UgJyArIHM7XG4gICAgfVxuXG4gICAgdmFyIG5zID0gdGhpcy5zdGFydCwgbmUgPSB0aGlzLmVuZDtcbiAgICBpZiAobnMgJiYgcykge1xuICAgICAgICBucyA9IG5zICsgcztcbiAgICB9IGVsc2Uge1xuICAgICAgICBucyA9IHMgfHwgbnM7XG4gICAgfVxuICAgIGlmIChsICYmIG5zKSB7XG4gICAgICAgIG5lID0gbnMgKyBsIC0gMTtcbiAgICB9IGVsc2Uge1xuICAgICAgICBuZSA9IG5lIHx8IGwgLSAxO1xuICAgIH1cbiAgICByZXR1cm4gbmV3IFVSTEZldGNoYWJsZSh0aGlzLnVybCwgbnMsIG5lLCB0aGlzLm9wdHMpO1xufVxuXG52YXIgc2VlZD0wO1xudmFyIGlzU2FmYXJpID0gbmF2aWdhdG9yLnVzZXJBZ2VudC5pbmRleE9mKCdTYWZhcmknKSA+PSAwICYmIG5hdmlnYXRvci51c2VyQWdlbnQuaW5kZXhPZignQ2hyb21lJykgPCAwIDtcblxuVVJMRmV0Y2hhYmxlLnByb3RvdHlwZS5mZXRjaEFzVGV4dCA9IGZ1bmN0aW9uKGNhbGxiYWNrKSB7XG4gICAgdmFyIHJlcSA9IG5ldyBYTUxIdHRwUmVxdWVzdCgpO1xuICAgIHZhciBsZW5ndGg7XG4gICAgdmFyIHVybCA9IHRoaXMudXJsO1xuICAgIGlmIChpc1NhZmFyaSB8fCB0aGlzLm9wdHMuc2FsdCkge1xuICAgICAgICB1cmwgPSB1cmwgKyAnP3NhbHQ9JyArIGI2NF9zaGExKCcnICsgRGF0ZS5ub3coKSArICcsJyArICgrK3NlZWQpKTtcbiAgICB9XG4gICAgcmVxLm9wZW4oJ0dFVCcsIHVybCwgdHJ1ZSk7XG5cbiAgICBpZiAodGhpcy5lbmQpIHtcbiAgICAgICAgaWYgKHRoaXMuZW5kIC0gdGhpcy5zdGFydCA+IDEwMDAwMDAwMCkge1xuICAgICAgICAgICAgdGhyb3cgJ01vbnN0ZXIgZmV0Y2ghJztcbiAgICAgICAgfVxuICAgICAgICByZXEuc2V0UmVxdWVzdEhlYWRlcignUmFuZ2UnLCAnYnl0ZXM9JyArIHRoaXMuc3RhcnQgKyAnLScgKyB0aGlzLmVuZCk7XG4gICAgICAgIGxlbmd0aCA9IHRoaXMuZW5kIC0gdGhpcy5zdGFydCArIDE7XG4gICAgfVxuXG4gICAgcmVxLm9ucmVhZHlzdGF0ZWNoYW5nZSA9IGZ1bmN0aW9uKCkge1xuICAgICAgICBpZiAocmVxLnJlYWR5U3RhdGUgPT0gNCkge1xuICAgICAgICAgICAgaWYgKHJlcS5zdGF0dXMgPT0gMjAwIHx8IHJlcS5zdGF0dXMgPT0gMjA2KSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKHJlcS5yZXNwb25zZVRleHQpO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9O1xuICAgIGlmICh0aGlzLm9wdHMuY3JlZGVudGlhbHMpIHtcbiAgICAgICAgcmVxLndpdGhDcmVkZW50aWFscyA9IHRydWU7XG4gICAgfVxuICAgIHJlcS5zZW5kKCcnKTtcbn1cblxuVVJMRmV0Y2hhYmxlLnByb3RvdHlwZS5zYWx0ZWQgPSBmdW5jdGlvbigpIHtcbiAgICB2YXIgbyA9IHNoYWxsb3dDb3B5KHRoaXMub3B0cyk7XG4gICAgby5zYWx0ID0gdHJ1ZTtcbiAgICByZXR1cm4gbmV3IFVSTEZldGNoYWJsZSh0aGlzLnVybCwgdGhpcy5zdGFydCwgdGhpcy5lbmQsIG8pO1xufVxuXG5VUkxGZXRjaGFibGUucHJvdG90eXBlLmZldGNoID0gZnVuY3Rpb24oY2FsbGJhY2ssIGF0dGVtcHQsIHRydW5jYXRlZExlbmd0aCkge1xuICAgIHZhciB0aGlzQiA9IHRoaXM7XG5cbiAgICBhdHRlbXB0ID0gYXR0ZW1wdCB8fCAxO1xuICAgIGlmIChhdHRlbXB0ID4gMykge1xuICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgfVxuXG4gICAgdmFyIHJlcSA9IG5ldyBYTUxIdHRwUmVxdWVzdCgpO1xuICAgIHZhciBsZW5ndGg7XG4gICAgdmFyIHVybCA9IHRoaXMudXJsO1xuICAgIGlmIChpc1NhZmFyaSB8fCB0aGlzLm9wdHMuc2FsdCkge1xuICAgICAgICB1cmwgPSB1cmwgKyAnP3NhbHQ9JyArIGI2NF9zaGExKCcnICsgRGF0ZS5ub3coKSArICcsJyArICgrK3NlZWQpKTtcbiAgICB9XG4gICAgcmVxLm9wZW4oJ0dFVCcsIHVybCwgdHJ1ZSk7XG4gICAgcmVxLm92ZXJyaWRlTWltZVR5cGUoJ3RleHQvcGxhaW47IGNoYXJzZXQ9eC11c2VyLWRlZmluZWQnKTtcbiAgICBpZiAodGhpcy5lbmQpIHtcbiAgICAgICAgaWYgKHRoaXMuZW5kIC0gdGhpcy5zdGFydCA+IDEwMDAwMDAwMCkge1xuICAgICAgICAgICAgdGhyb3cgJ01vbnN0ZXIgZmV0Y2ghJztcbiAgICAgICAgfVxuICAgICAgICByZXEuc2V0UmVxdWVzdEhlYWRlcignUmFuZ2UnLCAnYnl0ZXM9JyArIHRoaXMuc3RhcnQgKyAnLScgKyB0aGlzLmVuZCk7XG4gICAgICAgIGxlbmd0aCA9IHRoaXMuZW5kIC0gdGhpcy5zdGFydCArIDE7XG4gICAgfVxuICAgIHJlcS5yZXNwb25zZVR5cGUgPSAnYXJyYXlidWZmZXInO1xuICAgIHJlcS5vbnJlYWR5c3RhdGVjaGFuZ2UgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgaWYgKHJlcS5yZWFkeVN0YXRlID09IDQpIHtcbiAgICAgICAgICAgIGlmIChyZXEuc3RhdHVzID09IDIwMCB8fCByZXEuc3RhdHVzID09IDIwNikge1xuICAgICAgICAgICAgICAgIGlmIChyZXEucmVzcG9uc2UpIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGJsID0gcmVxLnJlc3BvbnNlLmJ5dGVMZW5ndGg7XG4gICAgICAgICAgICAgICAgICAgIGlmIChsZW5ndGggJiYgbGVuZ3RoICE9IGJsICYmICghdHJ1bmNhdGVkTGVuZ3RoIHx8IGJsICE9IHRydW5jYXRlZExlbmd0aCkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiB0aGlzQi5mZXRjaChjYWxsYmFjaywgYXR0ZW1wdCArIDEsIGJsKTtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhyZXEucmVzcG9uc2UpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfSBlbHNlIGlmIChyZXEubW96UmVzcG9uc2VBcnJheUJ1ZmZlcikge1xuICAgICAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2socmVxLm1velJlc3BvbnNlQXJyYXlCdWZmZXIpO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciByID0gcmVxLnJlc3BvbnNlVGV4dDtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGxlbmd0aCAmJiBsZW5ndGggIT0gci5sZW5ndGggJiYgKCF0cnVuY2F0ZWRMZW5ndGggfHwgci5sZW5ndGggIT0gdHJ1bmNhdGVkTGVuZ3RoKSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIHRoaXNCLmZldGNoKGNhbGxiYWNrLCBhdHRlbXB0ICsgMSwgci5sZW5ndGgpO1xuICAgICAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGJzdHJpbmdUb0J1ZmZlcihyZXEucmVzcG9uc2VUZXh0KSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIHJldHVybiB0aGlzQi5mZXRjaChjYWxsYmFjaywgYXR0ZW1wdCArIDEpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfTtcbiAgICBpZiAodGhpcy5vcHRzLmNyZWRlbnRpYWxzKSB7XG4gICAgICAgIHJlcS53aXRoQ3JlZGVudGlhbHMgPSB0cnVlO1xuICAgIH1cbiAgICByZXEuc2VuZCgnJyk7XG59XG5cbmZ1bmN0aW9uIGJzdHJpbmdUb0J1ZmZlcihyZXN1bHQpIHtcbiAgICBpZiAoIXJlc3VsdCkge1xuICAgICAgICByZXR1cm4gbnVsbDtcbiAgICB9XG5cbiAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShyZXN1bHQubGVuZ3RoKTtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGJhLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIGJhW2ldID0gcmVzdWx0LmNoYXJDb2RlQXQoaSk7XG4gICAgfVxuICAgIHJldHVybiBiYS5idWZmZXI7XG59XG5cbi8vIFJlYWQgZnJvbSBVaW50OEFycmF5XG5cbihmdW5jdGlvbihnbG9iYWwpIHtcbiAgICB2YXIgY29udmVydEJ1ZmZlciA9IG5ldyBBcnJheUJ1ZmZlcig4KTtcbiAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShjb252ZXJ0QnVmZmVyKTtcbiAgICB2YXIgZmEgPSBuZXcgRmxvYXQzMkFycmF5KGNvbnZlcnRCdWZmZXIpO1xuXG5cbiAgICBnbG9iYWwucmVhZEZsb2F0ID0gZnVuY3Rpb24oYnVmLCBvZmZzZXQpIHtcbiAgICAgICAgYmFbMF0gPSBidWZbb2Zmc2V0XTtcbiAgICAgICAgYmFbMV0gPSBidWZbb2Zmc2V0KzFdO1xuICAgICAgICBiYVsyXSA9IGJ1ZltvZmZzZXQrMl07XG4gICAgICAgIGJhWzNdID0gYnVmW29mZnNldCszXTtcbiAgICAgICAgcmV0dXJuIGZhWzBdO1xuICAgIH07XG4gfSh0aGlzKSk7XG5cbmZ1bmN0aW9uIHJlYWRJbnQ2NChiYSwgb2Zmc2V0KSB7XG4gICAgcmV0dXJuIChiYVtvZmZzZXQgKyA3XSA8PCAyNCkgfCAoYmFbb2Zmc2V0ICsgNl0gPDwgMTYpIHwgKGJhW29mZnNldCArIDVdIDw8IDgpIHwgKGJhW29mZnNldCArIDRdKTtcbn1cblxuZnVuY3Rpb24gcmVhZEludChiYSwgb2Zmc2V0KSB7XG4gICAgcmV0dXJuIChiYVtvZmZzZXQgKyAzXSA8PCAyNCkgfCAoYmFbb2Zmc2V0ICsgMl0gPDwgMTYpIHwgKGJhW29mZnNldCArIDFdIDw8IDgpIHwgKGJhW29mZnNldF0pO1xufVxuXG5mdW5jdGlvbiByZWFkU2hvcnQoYmEsIG9mZnNldCkge1xuICAgIHJldHVybiAoYmFbb2Zmc2V0ICsgMV0gPDwgOCkgfCAoYmFbb2Zmc2V0XSk7XG59XG5cbmZ1bmN0aW9uIHJlYWRCeXRlKGJhLCBvZmZzZXQpIHtcbiAgICByZXR1cm4gYmFbb2Zmc2V0XTtcbn1cblxuZnVuY3Rpb24gcmVhZEludEJFKGJhLCBvZmZzZXQpIHtcbiAgICByZXR1cm4gKGJhW29mZnNldF0gPDwgMjQpIHwgKGJhW29mZnNldCArIDFdIDw8IDE2KSB8IChiYVtvZmZzZXQgKyAyXSA8PCA4KSB8IChiYVtvZmZzZXQgKyAzXSk7XG59XG5cbi8vIEV4cG9ydHMgaWYgd2UgYXJlIGJlaW5nIHVzZWQgYXMgYSBtb2R1bGVcblxuaWYgKHR5cGVvZihtb2R1bGUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIG1vZHVsZS5leHBvcnRzID0ge1xuICAgICAgICBCbG9iRmV0Y2hhYmxlOiBCbG9iRmV0Y2hhYmxlLFxuICAgICAgICBVUkxGZXRjaGFibGU6IFVSTEZldGNoYWJsZSxcblxuICAgICAgICByZWFkSW50OiByZWFkSW50LFxuICAgICAgICByZWFkSW50QkU6IHJlYWRJbnRCRSxcbiAgICAgICAgcmVhZEludDY0OiByZWFkSW50NjQsXG4gICAgICAgIHJlYWRTaG9ydDogcmVhZFNob3J0LFxuICAgICAgICByZWFkQnl0ZTogcmVhZEJ5dGUsXG4gICAgICAgIHJlYWRGbG9hdDogdGhpcy5yZWFkRmxvYXRcbiAgICB9XG59XG4iLCIvKiAtKi0gbW9kZTogamF2YXNjcmlwdDsgYy1iYXNpYy1vZmZzZXQ6IDQ7IGluZGVudC10YWJzLW1vZGU6IG5pbCAtKi0gKi9cblxuLy8gXG4vLyBEYWxsaWFuY2UgR2Vub21lIEV4cGxvcmVyXG4vLyAoYykgVGhvbWFzIERvd24gMjAwNi0yMDEwXG4vL1xuLy8gY29sb3IuanNcbi8vXG5cblwidXNlIHN0cmljdFwiO1xuXG5mdW5jdGlvbiBEQ29sb3VyKHJlZCwgZ3JlZW4sIGJsdWUsIG5hbWUpIHtcbiAgICB0aGlzLnJlZCA9IHJlZHwwO1xuICAgIHRoaXMuZ3JlZW4gPSBncmVlbnwwO1xuICAgIHRoaXMuYmx1ZSA9IGJsdWV8MDtcbiAgICBpZiAobmFtZSkge1xuICAgICAgICB0aGlzLm5hbWUgPSBuYW1lO1xuICAgIH1cbn1cblxuRENvbG91ci5wcm90b3R5cGUudG9TdmdTdHJpbmcgPSBmdW5jdGlvbigpIHtcbiAgICBpZiAoIXRoaXMubmFtZSkge1xuICAgICAgICB0aGlzLm5hbWUgPSBcInJnYihcIiArIHRoaXMucmVkICsgXCIsXCIgKyB0aGlzLmdyZWVuICsgXCIsXCIgKyB0aGlzLmJsdWUgKyBcIilcIjtcbiAgICB9XG5cbiAgICByZXR1cm4gdGhpcy5uYW1lO1xufVxuXG5mdW5jdGlvbiBoZXgyKHgpIHtcbiAgICB2YXIgeSA9ICcwMCcgKyB4LnRvU3RyaW5nKDE2KTtcbiAgICByZXR1cm4geS5zdWJzdHJpbmcoeS5sZW5ndGggLSAyKTtcbn1cblxuRENvbG91ci5wcm90b3R5cGUudG9IZXhTdHJpbmcgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gJyMnICsgaGV4Mih0aGlzLnJlZCkgKyBoZXgyKHRoaXMuZ3JlZW4pICsgaGV4Mih0aGlzLmJsdWUpO1xufVxuXG52YXIgcGFsZXR0ZSA9IHtcbiAgICByZWQ6IG5ldyBEQ29sb3VyKDI1NSwgMCwgMCwgJ3JlZCcpLFxuICAgIGdyZWVuOiBuZXcgRENvbG91cigwLCAyNTUsIDAsICdncmVlbicpLFxuICAgIGJsdWU6IG5ldyBEQ29sb3VyKDAsIDAsIDI1NSwgJ2JsdWUnKSxcbiAgICB5ZWxsb3c6IG5ldyBEQ29sb3VyKDI1NSwgMjU1LCAwLCAneWVsbG93JyksXG4gICAgd2hpdGU6IG5ldyBEQ29sb3VyKDI1NSwgMjU1LCAyNTUsICd3aGl0ZScpLFxuICAgIGJsYWNrOiBuZXcgRENvbG91cigwLCAwLCAwLCAnYmxhY2snKSxcbiAgICBncmF5OiBuZXcgRENvbG91cigxODAsIDE4MCwgMTgwLCAnZ3JheScpLFxuICAgIGdyZXk6IG5ldyBEQ29sb3VyKDE4MCwgMTgwLCAxODAsICdncmV5JyksXG4gICAgbGlnaHRza3libHVlOiBuZXcgRENvbG91cigxMzUsIDIwNiwgMjUwLCAnbGlnaHRza3libHVlJyksXG4gICAgbGlnaHRzYWxtb246IG5ldyBEQ29sb3VyKDI1NSwgMTYwLCAxMjIsICdsaWdodHNhbG1vbicpLFxuICAgIGhvdHBpbms6IG5ldyBEQ29sb3VyKDI1NSwgMTA1LCAxODAsICdob3RwaW5rJylcbn07XG5cbnZhciBDT0xPUl9SRSA9IG5ldyBSZWdFeHAoJ14jKFswLTlBLUZhLWZdezJ9KShbMC05QS1GYS1mXXsyfSkoWzAtOUEtRmEtZl17Mn0pJCcpO1xudmFyIENTU19DT0xPUl9SRSA9IC9yZ2JcXCgoWzAtOV0rKSwoWzAtOV0rKSwoWzAtOV0rKVxcKS9cblxuZnVuY3Rpb24gZGFzQ29sb3VyRm9yTmFtZShuYW1lKSB7XG4gICAgdmFyIGMgPSBwYWxldHRlW25hbWVdO1xuICAgIGlmICghYykge1xuICAgICAgICB2YXIgbWF0Y2ggPSBDT0xPUl9SRS5leGVjKG5hbWUpO1xuICAgICAgICBpZiAobWF0Y2gpIHtcbiAgICAgICAgICAgIGMgPSBuZXcgRENvbG91cigoJzB4JyArIG1hdGNoWzFdKXwwLCAoJzB4JyArIG1hdGNoWzJdKXwwLCAoJzB4JyArIG1hdGNoWzNdKXwwLCBuYW1lKTtcbiAgICAgICAgICAgIHBhbGV0dGVbbmFtZV0gPSBjO1xuICAgICAgICB9IGVsc2Uge1xuICAgIFx0ICAgIG1hdGNoID0gQ1NTX0NPTE9SX1JFLmV4ZWMobmFtZSk7XG4gICAgXHQgICAgaWYgKG1hdGNoKSB7XG4gICAgICAgIFx0XHRjID0gbmV3IERDb2xvdXIobWF0Y2hbMV18MCwgbWF0Y2hbMl18MCwgbWF0Y2hbM118MCwgbmFtZSk7XG4gICAgICAgIFx0XHRwYWxldHRlW25hbWVdID0gYztcblx0ICAgICAgIH0gZWxzZSB7XG5cdFx0ICAgICAgY29uc29sZS5sb2coXCJjb3VsZG4ndCBoYW5kbGUgY29sb3I6IFwiICsgbmFtZSk7XG5cdFx0ICAgICAgYyA9IHBhbGV0dGUuYmxhY2s7XG5cdFx0ICAgICAgcGFsZXR0ZVtuYW1lXSA9IGM7XG5cdCAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuIGM7XG59XG5cbmZ1bmN0aW9uIG1ha2VDb2xvdXJTdGVwcyhzdGVwcywgc3RvcHMsIGNvbG91cnMpIHtcbiAgICB2YXIgZGNvbG91cnMgPSBbXTtcbiAgICBmb3IgKHZhciBjaSA9IDA7IGNpIDwgY29sb3Vycy5sZW5ndGg7ICsrY2kpIHtcbiAgICAgICAgZGNvbG91cnMucHVzaChkYXNDb2xvdXJGb3JOYW1lKGNvbG91cnNbY2ldKSk7XG4gICAgfVxuXG4gICAgdmFyIGdyYWQgPSBbXTtcbiAgU1RFUF9MT09QOlxuICAgIGZvciAodmFyIHNpID0gMDsgc2kgPCBzdGVwczsgKytzaSkge1xuICAgICAgICB2YXIgcnMgPSAoMS4wICogc2kpIC8gKHN0ZXBzLTEpO1xuICAgICAgICB2YXIgc2NvcmUgPSBzdG9wc1swXSArIChzdG9wc1tzdG9wcy5sZW5ndGggLTFdIC0gc3RvcHNbMF0pICogcnM7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgc3RvcHMubGVuZ3RoIC0gMTsgKytpKSB7XG4gICAgICAgICAgICBpZiAoc2NvcmUgPj0gc3RvcHNbaV0gJiYgc2NvcmUgPD0gc3RvcHNbaSsxXSkge1xuICAgICAgICAgICAgICAgIHZhciBmcmFjID0gKHNjb3JlIC0gc3RvcHNbaV0pIC8gKHN0b3BzW2krMV0gLSBzdG9wc1tpXSk7XG4gICAgICAgICAgICAgICAgdmFyIGNhID0gZGNvbG91cnNbaV07XG4gICAgICAgICAgICAgICAgdmFyIGNiID0gZGNvbG91cnNbaSsxXTtcblxuICAgICAgICAgICAgICAgIHZhciBmaWxsID0gbmV3IERDb2xvdXIoXG4gICAgICAgICAgICAgICAgICAgICgoY2EucmVkICogKDEuMCAtIGZyYWMpKSArIChjYi5yZWQgKiBmcmFjKSl8MCxcbiAgICAgICAgICAgICAgICAgICAgKChjYS5ncmVlbiAqICgxLjAgLSBmcmFjKSkgKyAoY2IuZ3JlZW4gKiBmcmFjKSl8MCxcbiAgICAgICAgICAgICAgICAgICAgKChjYS5ibHVlICogKDEuMCAtIGZyYWMpKSArIChjYi5ibHVlICogZnJhYykpfDBcbiAgICAgICAgICAgICAgICApLnRvU3ZnU3RyaW5nKCk7XG4gICAgICAgICAgICAgICAgZ3JhZC5wdXNoKGZpbGwpO1xuXG4gICAgICAgICAgICAgICAgY29udGludWUgU1RFUF9MT09QO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHRocm93ICdCYWQgc3RlcCc7XG4gICAgfVxuXG4gICAgcmV0dXJuIGdyYWQ7XG59XG5cbmZ1bmN0aW9uIG1ha2VHcmFkaWVudChzdGVwcywgY29sb3IxLCBjb2xvcjIsIGNvbG9yMykge1xuICAgIGlmIChjb2xvcjMpIHtcbiAgICAgICAgcmV0dXJuIG1ha2VDb2xvdXJTdGVwcyhzdGVwcywgWzAsIDAuNSwgMV0sIFtjb2xvcjEsIGNvbG9yMiwgY29sb3IzXSk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgcmV0dXJuIG1ha2VDb2xvdXJTdGVwcyhzdGVwcywgWzAsIDFdLCBbY29sb3IxLCBjb2xvcjJdKTtcbiAgICB9XG59XG5cbmlmICh0eXBlb2YobW9kdWxlKSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgICBtb2R1bGUuZXhwb3J0cyA9IHtcbiAgICAgICAgbWFrZUNvbG91clN0ZXBzOiBtYWtlQ29sb3VyU3RlcHMsXG4gICAgICAgIG1ha2VHcmFkaWVudDogbWFrZUdyYWRpZW50LFxuICAgICAgICBkYXNDb2xvdXJGb3JOYW1lOiBkYXNDb2xvdXJGb3JOYW1lXG4gICAgfTtcbn1cbiIsIi8qIC0qLSBtb2RlOiBqYXZhc2NyaXB0OyBjLWJhc2ljLW9mZnNldDogNDsgaW5kZW50LXRhYnMtbW9kZTogbmlsIC0qLSAqL1xuXG4vLyBcbi8vIERhbGxpYW5jZSBHZW5vbWUgRXhwbG9yZXJcbi8vIChjKSBUaG9tYXMgRG93biAyMDA2LTIwMTBcbi8vXG4vLyBkYXMuanM6IHF1ZXJpZXMgYW5kIGxvdy1sZXZlbCBkYXRhIG1vZGVsLlxuLy9cblxuXCJ1c2Ugc3RyaWN0XCI7XG5cbmlmICh0eXBlb2YocmVxdWlyZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgdmFyIHV0aWxzID0gcmVxdWlyZSgnLi91dGlscycpO1xuICAgIHZhciBzaGFsbG93Q29weSA9IHV0aWxzLnNoYWxsb3dDb3B5O1xuICAgIHZhciBwdXNobyA9IHV0aWxzLnB1c2hvO1xuXG4gICAgdmFyIGNvbG9yID0gcmVxdWlyZSgnLi9jb2xvcicpO1xuICAgIHZhciBtYWtlQ29sb3VyU3RlcHMgPSBjb2xvci5tYWtlQ29sb3VyU3RlcHM7XG59XG5cbnZhciBkYXNMaWJFcnJvckhhbmRsZXIgPSBmdW5jdGlvbihlcnJNc2cpIHtcbiAgICBhbGVydChlcnJNc2cpO1xufVxudmFyIGRhc0xpYlJlcXVlc3RRdWV1ZSA9IG5ldyBBcnJheSgpO1xuXG5cblxuZnVuY3Rpb24gREFTU2VnbWVudChuYW1lLCBzdGFydCwgZW5kLCBkZXNjcmlwdGlvbikge1xuICAgIHRoaXMubmFtZSA9IG5hbWU7XG4gICAgdGhpcy5zdGFydCA9IHN0YXJ0O1xuICAgIHRoaXMuZW5kID0gZW5kO1xuICAgIHRoaXMuZGVzY3JpcHRpb24gPSBkZXNjcmlwdGlvbjtcbn1cbkRBU1NlZ21lbnQucHJvdG90eXBlLnRvU3RyaW5nID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMubmFtZSArICc6JyArIHRoaXMuc3RhcnQgKyAnLi4nICsgdGhpcy5lbmQ7XG59O1xuREFTU2VnbWVudC5wcm90b3R5cGUuaXNCb3VuZGVkID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuc3RhcnQgJiYgdGhpcy5lbmQ7XG59XG5EQVNTZWdtZW50LnByb3RvdHlwZS50b0RBU1F1ZXJ5ID0gZnVuY3Rpb24oKSB7XG4gICAgdmFyIHEgPSAnc2VnbWVudD0nICsgdGhpcy5uYW1lO1xuICAgIGlmICh0aGlzLnN0YXJ0ICYmIHRoaXMuZW5kKSB7XG4gICAgICAgIHEgKz0gKCc6JyArIHRoaXMuc3RhcnQgKyAnLCcgKyB0aGlzLmVuZCk7XG4gICAgfVxuICAgIHJldHVybiBxO1xufVxuXG5cbmZ1bmN0aW9uIERBU1NvdXJjZShhMSwgYTIpIHtcbiAgICB2YXIgb3B0aW9ucztcbiAgICBpZiAodHlwZW9mIGExID09ICdzdHJpbmcnKSB7XG4gICAgICAgIHRoaXMudXJpID0gYTE7XG4gICAgICAgIG9wdGlvbnMgPSBhMiB8fCB7fTtcbiAgICB9IGVsc2Uge1xuICAgICAgICBvcHRpb25zID0gYTEgfHwge307XG4gICAgfVxuICAgIGZvciAodmFyIGsgaW4gb3B0aW9ucykge1xuICAgICAgICBpZiAodHlwZW9mKG9wdGlvbnNba10pICE9ICdmdW5jdGlvbicpIHtcbiAgICAgICAgICAgIHRoaXNba10gPSBvcHRpb25zW2tdO1xuICAgICAgICB9XG4gICAgfVxuXG5cbiAgICBpZiAoIXRoaXMuY29vcmRzKSB7XG4gICAgICAgIHRoaXMuY29vcmRzID0gW107XG4gICAgfVxuICAgIGlmICghdGhpcy5wcm9wcykge1xuICAgICAgICB0aGlzLnByb3BzID0ge307XG4gICAgfVxuXG4gICAgdGhpcy5kYXNCYXNlVVJJID0gdGhpcy51cmk7XG4gICAgaWYgKHRoaXMuZGFzQmFzZVVSSSAmJiB0aGlzLmRhc0Jhc2VVUkkuc3Vic3RyKHRoaXMudXJpLmxlbmd0aCAtIDEpICE9ICcvJykge1xuICAgICAgICB0aGlzLmRhc0Jhc2VVUkkgPSB0aGlzLmRhc0Jhc2VVUkkgKyAnLyc7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBEQVNDb29yZHMoKSB7XG59XG5cbmZ1bmN0aW9uIGNvb3Jkc01hdGNoKGMxLCBjMikge1xuICAgIHJldHVybiBjMS50YXhvbiA9PSBjMi50YXhvbiAmJiBjMS5hdXRoID09IGMyLmF1dGggJiYgYzEudmVyc2lvbiA9PSBjMi52ZXJzaW9uO1xufVxuXG4vL1xuLy8gREFTIDEuNiBlbnRyeV9wb2ludHMgY29tbWFuZFxuLy9cblxuREFTU291cmNlLnByb3RvdHlwZS5lbnRyeVBvaW50cyA9IGZ1bmN0aW9uKGNhbGxiYWNrKSB7XG4gICAgdmFyIGRhc1VSSSA9IHRoaXMuZGFzQmFzZVVSSSArICdlbnRyeV9wb2ludHMnO1xuICAgIHRoaXMuZG9Dcm9zc0RvbWFpblJlcXVlc3QoZGFzVVJJLCBmdW5jdGlvbihyZXNwb25zZVhNTCkge1xuICAgICAgICAgICAgaWYgKCFyZXNwb25zZVhNTCkge1xuICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhbXSk7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICB2YXIgZW50cnlQb2ludHMgPSBuZXcgQXJyYXkoKTtcbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICB2YXIgc2VncyA9IHJlc3BvbnNlWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdTRUdNRU5UJyk7XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBzZWdzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWcgPSBzZWdzW2ldO1xuICAgICAgICAgICAgICAgICAgICB2YXIgc2VnSWQgPSBzZWcuZ2V0QXR0cmlidXRlKCdpZCcpO1xuICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgdmFyIHNlZ1NpemUgPSBzZWcuZ2V0QXR0cmlidXRlKCdzaXplJyk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWdNaW4sIHNlZ01heDtcbiAgICAgICAgICAgICAgICAgICAgaWYgKHNlZ1NpemUpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHNlZ01pbiA9IDE7IHNlZ01heCA9IHNlZ1NpemV8MDtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHNlZ01pbiA9IHNlZy5nZXRBdHRyaWJ1dGUoJ3N0YXJ0Jyk7XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoc2VnTWluKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgc2VnTWluIHw9IDA7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICBzZWdNYXggPSBzZWcuZ2V0QXR0cmlidXRlKCdzdG9wJyk7XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoc2VnTWF4KSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgc2VnTWF4IHw9IDA7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgdmFyIHNlZ0Rlc2MgPSBudWxsO1xuICAgICAgICAgICAgICAgICAgICBpZiAoc2VnLmZpcnN0Q2hpbGQpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHNlZ0Rlc2MgPSBzZWcuZmlyc3RDaGlsZC5ub2RlVmFsdWU7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgZW50cnlQb2ludHMucHVzaChuZXcgREFTU2VnbWVudChzZWdJZCwgc2VnTWluLCBzZWdNYXgsIHNlZ0Rlc2MpKTtcbiAgICAgICAgICAgICAgICB9ICAgICAgICAgIFxuICAgICAgICAgICAgICAgY2FsbGJhY2soZW50cnlQb2ludHMpO1xuICAgIH0pOyAgICAgICAgIFxufVxuXG4vL1xuLy8gREFTIDEuNiBzZXF1ZW5jZSBjb21tYW5kXG4vLyBEbyB3ZSBuZWVkIGFuIG9wdGlvbiB0byBmYWxsIGJhY2sgdG8gdGhlIGRuYSBjb21tYW5kP1xuLy9cblxuZnVuY3Rpb24gREFTU2VxdWVuY2UobmFtZSwgc3RhcnQsIGVuZCwgYWxwaGEsIHNlcSkge1xuICAgIHRoaXMubmFtZSA9IG5hbWU7XG4gICAgdGhpcy5zdGFydCA9IHN0YXJ0O1xuICAgIHRoaXMuZW5kID0gZW5kO1xuICAgIHRoaXMuYWxwaGFiZXQgPSBhbHBoYTtcbiAgICB0aGlzLnNlcSA9IHNlcTtcbn1cblxuREFTU291cmNlLnByb3RvdHlwZS5zZXF1ZW5jZSA9IGZ1bmN0aW9uKHNlZ21lbnQsIGNhbGxiYWNrKSB7XG4gICAgdmFyIGRhc1VSSSA9IHRoaXMuZGFzQmFzZVVSSSArICdzZXF1ZW5jZT8nICsgc2VnbWVudC50b0RBU1F1ZXJ5KCk7XG4gICAgdGhpcy5kb0Nyb3NzRG9tYWluUmVxdWVzdChkYXNVUkksIGZ1bmN0aW9uKHJlc3BvbnNlWE1MKSB7XG4gICAgICAgIGlmICghcmVzcG9uc2VYTUwpIHtcbiAgICAgICAgICAgIGNhbGxiYWNrKFtdKTtcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICB2YXIgc2VxcyA9IG5ldyBBcnJheSgpO1xuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIHZhciBzZWdzID0gcmVzcG9uc2VYTUwuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1NFUVVFTkNFJyk7XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBzZWdzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWcgPSBzZWdzW2ldO1xuICAgICAgICAgICAgICAgICAgICB2YXIgc2VnSWQgPSBzZWcuZ2V0QXR0cmlidXRlKCdpZCcpO1xuICAgICAgICAgICAgICAgICAgICB2YXIgc2VnTWluID0gc2VnLmdldEF0dHJpYnV0ZSgnc3RhcnQnKTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHNlZ01heCA9IHNlZy5nZXRBdHRyaWJ1dGUoJ3N0b3AnKTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHNlZ0FscGhhID0gJ0ROQSc7XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWdTZXEgPSBudWxsO1xuICAgICAgICAgICAgICAgICAgICBpZiAoc2VnLmZpcnN0Q2hpbGQpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciByYXdTZXEgPSBzZWcuZmlyc3RDaGlsZC5ub2RlVmFsdWU7XG4gICAgICAgICAgICAgICAgICAgICAgICBzZWdTZXEgPSAnJztcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBpZHggPSAwO1xuICAgICAgICAgICAgICAgICAgICAgICAgd2hpbGUgKHRydWUpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgc3BhY2UgPSByYXdTZXEuaW5kZXhPZignXFxuJywgaWR4KTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoc3BhY2UgPj0gMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBzZWdTZXEgKz0gcmF3U2VxLnN1YnN0cmluZyhpZHgsIHNwYWNlKS50b1VwcGVyQ2FzZSgpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZHggPSBzcGFjZSArIDE7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgc2VnU2VxICs9IHJhd1NlcS5zdWJzdHJpbmcoaWR4KS50b1VwcGVyQ2FzZSgpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgc2Vxcy5wdXNoKG5ldyBEQVNTZXF1ZW5jZShzZWdJZCwgc2VnTWluLCBzZWdNYXgsIHNlZ0FscGhhLCBzZWdTZXEpKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgY2FsbGJhY2soc2Vxcyk7XG4gICAgICAgIH1cbiAgICB9KTtcbn1cblxuLy9cbi8vIERBUyAxLjYgZmVhdHVyZXMgY29tbWFuZFxuLy9cblxuZnVuY3Rpb24gREFTRmVhdHVyZSgpIHtcbn1cblxuZnVuY3Rpb24gREFTR3JvdXAoaWQpIHtcbiAgICBpZiAoaWQpXG4gICAgICAgIHRoaXMuaWQgPSBpZDtcbn1cblxuZnVuY3Rpb24gREFTTGluayhkZXNjLCB1cmkpIHtcbiAgICB0aGlzLmRlc2MgPSBkZXNjO1xuICAgIHRoaXMudXJpID0gdXJpO1xufVxuXG5EQVNTb3VyY2UucHJvdG90eXBlLmZlYXR1cmVzID0gZnVuY3Rpb24oc2VnbWVudCwgb3B0aW9ucywgY2FsbGJhY2spIHtcbiAgICBvcHRpb25zID0gb3B0aW9ucyB8fCB7fTtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuXG4gICAgdmFyIGRhc1VSSTtcbiAgICBpZiAodGhpcy5mZWF0dXJlc191cmkpIHtcbiAgICAgICAgZGFzVVJJID0gdGhpcy5mZWF0dXJlc191cmk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdmFyIGZpbHRlcnMgPSBbXTtcblxuICAgICAgICBpZiAoc2VnbWVudCkge1xuICAgICAgICAgICAgZmlsdGVycy5wdXNoKHNlZ21lbnQudG9EQVNRdWVyeSgpKTtcbiAgICAgICAgfSBlbHNlIGlmIChvcHRpb25zLmdyb3VwKSB7XG4gICAgICAgICAgICB2YXIgZyA9IG9wdGlvbnMuZ3JvdXA7XG4gICAgICAgICAgICBpZiAodHlwZW9mIGcgPT0gJ3N0cmluZycpIHtcbiAgICAgICAgICAgICAgICBmaWx0ZXJzLnB1c2goJ2dyb3VwX2lkPScgKyBnKTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgZ2kgPSAwOyBnaSA8IGcubGVuZ3RoOyArK2dpKSB7XG4gICAgICAgICAgICAgICAgICAgIGZpbHRlcnMucHVzaCgnZ3JvdXBfaWQ9JyArIGdbZ2ldKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBpZiAob3B0aW9ucy5hZGphY2VudCkge1xuICAgICAgICAgICAgdmFyIGFkaiA9IG9wdGlvbnMuYWRqYWNlbnQ7XG4gICAgICAgICAgICBpZiAodHlwZW9mIGFkaiA9PSAnc3RyaW5nJykge1xuICAgICAgICAgICAgICAgIGFkaiA9IFthZGpdO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZm9yICh2YXIgYWkgPSAwOyBhaSA8IGFkai5sZW5ndGg7ICsrYWkpIHtcbiAgICAgICAgICAgICAgICBmaWx0ZXJzLnB1c2goJ2FkamFjZW50PScgKyBhZGpbYWldKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGlmIChvcHRpb25zLnR5cGUpIHtcbiAgICAgICAgICAgIGlmICh0eXBlb2Ygb3B0aW9ucy50eXBlID09ICdzdHJpbmcnKSB7XG4gICAgICAgICAgICAgICAgZmlsdGVycy5wdXNoKCd0eXBlPScgKyBvcHRpb25zLnR5cGUpO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICBmb3IgKHZhciB0aSA9IDA7IHRpIDwgb3B0aW9ucy50eXBlLmxlbmd0aDsgKyt0aSkge1xuICAgICAgICAgICAgICAgICAgICBmaWx0ZXJzLnB1c2goJ3R5cGU9JyArIG9wdGlvbnMudHlwZVt0aV0pO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBcbiAgICAgICAgaWYgKG9wdGlvbnMubWF4Ymlucykge1xuICAgICAgICAgICAgZmlsdGVycy5wdXNoKCdtYXhiaW5zPScgKyBvcHRpb25zLm1heGJpbnMpO1xuICAgICAgICB9XG4gICAgICAgIFxuICAgICAgICBpZiAoZmlsdGVycy5sZW5ndGggPiAwKSB7XG4gICAgICAgICAgICBkYXNVUkkgPSB0aGlzLmRhc0Jhc2VVUkkgKyAnZmVhdHVyZXM/JyArIGZpbHRlcnMuam9pbignOycpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgY2FsbGJhY2soW10sICdObyBmaWx0ZXJzIHNwZWNpZmllZCcpO1xuICAgICAgICB9XG4gICAgfSBcbiAgIFxuXG4gICAgdGhpcy5kb0Nyb3NzRG9tYWluUmVxdWVzdChkYXNVUkksIGZ1bmN0aW9uKHJlc3BvbnNlWE1MLCByZXEpIHtcbiAgICAgICAgaWYgKCFyZXNwb25zZVhNTCkge1xuICAgICAgICAgICAgdmFyIG1zZztcbiAgICAgICAgICAgIGlmIChyZXEuc3RhdHVzID09IDApIHtcbiAgICAgICAgICAgICAgICBtc2cgPSAnc2VydmVyIG1heSBub3Qgc3VwcG9ydCBDT1JTJztcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgbXNnID0gJ3N0YXR1cz0nICsgcmVxLnN0YXR1cztcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGNhbGxiYWNrKFtdLCAnRmFpbGVkIHJlcXVlc3Q6ICcgKyBtc2cpO1xuICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICB9XG4vKiAgICAgIGlmIChyZXEpIHtcbiAgICAgICAgICAgIHZhciBjYXBzID0gcmVxLmdldFJlc3BvbnNlSGVhZGVyKCdYLURBUy1DYXBhYmlsdGllcycpO1xuICAgICAgICAgICAgaWYgKGNhcHMpIHtcbiAgICAgICAgICAgICAgICBhbGVydChjYXBzKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSAqL1xuXG4gICAgICAgIHZhciBmZWF0dXJlcyA9IG5ldyBBcnJheSgpO1xuICAgICAgICB2YXIgc2VnbWVudE1hcCA9IHt9O1xuXG4gICAgICAgIHZhciBzZWdzID0gcmVzcG9uc2VYTUwuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1NFR01FTlQnKTtcbiAgICAgICAgZm9yICh2YXIgc2kgPSAwOyBzaSA8IHNlZ3MubGVuZ3RoOyArK3NpKSB7XG4gICAgICAgICAgICB2YXIgc2VnbWVudFhNTCA9IHNlZ3Nbc2ldO1xuICAgICAgICAgICAgdmFyIHNlZ21lbnRJRCA9IHNlZ21lbnRYTUwuZ2V0QXR0cmlidXRlKCdpZCcpO1xuICAgICAgICAgICAgc2VnbWVudE1hcFtzZWdtZW50SURdID0ge1xuICAgICAgICAgICAgICAgIG1pbjogc2VnbWVudFhNTC5nZXRBdHRyaWJ1dGUoJ3N0YXJ0JyksXG4gICAgICAgICAgICAgICAgbWF4OiBzZWdtZW50WE1MLmdldEF0dHJpYnV0ZSgnc3RvcCcpXG4gICAgICAgICAgICB9O1xuICAgICAgICAgICAgXG4gICAgICAgICAgICB2YXIgZmVhdHVyZVhNTHMgPSBzZWdtZW50WE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdGRUFUVVJFJyk7XG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGZlYXR1cmVYTUxzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICAgICAgdmFyIGZlYXR1cmUgPSBmZWF0dXJlWE1Mc1tpXTtcbiAgICAgICAgICAgICAgICB2YXIgZGFzRmVhdHVyZSA9IG5ldyBEQVNGZWF0dXJlKCk7XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5zZWdtZW50ID0gc2VnbWVudElEO1xuICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUuaWQgPSBmZWF0dXJlLmdldEF0dHJpYnV0ZSgnaWQnKTtcbiAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLmxhYmVsID0gZmVhdHVyZS5nZXRBdHRyaWJ1dGUoJ2xhYmVsJyk7XG5cblxuLypcbiAgICAgICAgICAgICAgICB2YXIgY2hpbGROb2RlcyA9IGZlYXR1cmUuY2hpbGROb2RlcztcbiAgICAgICAgICAgICAgICBmb3IgKHZhciBjID0gMDsgYyA8IGNoaWxkTm9kZXMubGVuZ3RoOyArK2MpIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGNuID0gY2hpbGROb2Rlc1tjXTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGNuLm5vZGVUeXBlID09IE5vZGUuRUxFTUVOVF9OT0RFKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIga2V5ID0gY24udGFnTmFtZTtcbiAgICAgICAgICAgICAgICAgICAgICAgIC8vdmFyIHZhbCA9IG51bGw7XG4gICAgICAgICAgICAgICAgICAgICAgICAvL2lmIChjbi5maXJzdENoaWxkKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAvLyAgIHZhbCA9IGNuLmZpcnN0Q2hpbGQubm9kZVZhbHVlO1xuICAgICAgICAgICAgICAgICAgICAgICAgLy99XG4gICAgICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlW2tleV0gPSAneCc7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9ICovXG5cblxuICAgICAgICAgICAgICAgIHZhciBzcG9zID0gZWxlbWVudFZhbHVlKGZlYXR1cmUsIFwiU1RBUlRcIik7XG4gICAgICAgICAgICAgICAgdmFyIGVwb3MgPSBlbGVtZW50VmFsdWUoZmVhdHVyZSwgXCJFTkRcIik7XG4gICAgICAgICAgICAgICAgaWYgKChzcG9zfDApID4gKGVwb3N8MCkpIHtcbiAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5taW4gPSBlcG9zfDA7XG4gICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUubWF4ID0gc3Bvc3wwO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUubWluID0gc3Bvc3wwO1xuICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLm1heCA9IGVwb3N8MDtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICB2YXIgdGVjID0gZmVhdHVyZS5nZXRFbGVtZW50c0J5VGFnTmFtZSgnVFlQRScpO1xuICAgICAgICAgICAgICAgICAgICBpZiAodGVjLmxlbmd0aCA+IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciB0ZSA9IHRlY1swXTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmICh0ZS5maXJzdENoaWxkKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS50eXBlID0gdGUuZmlyc3RDaGlsZC5ub2RlVmFsdWU7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLnR5cGVJZCA9IHRlLmdldEF0dHJpYnV0ZSgnaWQnKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUudHlwZUN2ID0gdGUuZ2V0QXR0cmlidXRlKCdjdklkJyk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZGFzRmVhdHVyZS50eXBlID0gZWxlbWVudFZhbHVlKGZlYXR1cmUsIFwiVFlQRVwiKTtcbiAgICAgICAgICAgICAgICBpZiAoIWRhc0ZlYXR1cmUudHlwZSAmJiBkYXNGZWF0dXJlLnR5cGVJZCkge1xuICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLnR5cGUgPSBkYXNGZWF0dXJlLnR5cGVJZDsgLy8gRklYTUU/XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUubWV0aG9kID0gZWxlbWVudFZhbHVlKGZlYXR1cmUsIFwiTUVUSE9EXCIpO1xuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIG9yaSA9IGVsZW1lbnRWYWx1ZShmZWF0dXJlLCBcIk9SSUVOVEFUSU9OXCIpO1xuICAgICAgICAgICAgICAgICAgICBpZiAoIW9yaSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgb3JpID0gJzAnO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUub3JpZW50YXRpb24gPSBvcmk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUuc2NvcmUgPSBlbGVtZW50VmFsdWUoZmVhdHVyZSwgXCJTQ09SRVwiKTtcbiAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLmxpbmtzID0gZGFzTGlua3NPZihmZWF0dXJlKTtcbiAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLm5vdGVzID0gZGFzTm90ZXNPZihmZWF0dXJlKTtcbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICB2YXIgZ3JvdXBzID0gZmVhdHVyZS5nZXRFbGVtZW50c0J5VGFnTmFtZShcIkdST1VQXCIpO1xuICAgICAgICAgICAgICAgIGZvciAodmFyIGdpICA9IDA7IGdpIDwgZ3JvdXBzLmxlbmd0aDsgKytnaSkge1xuICAgICAgICAgICAgICAgICAgICB2YXIgZ3JvdXBYTUwgPSBncm91cHNbZ2ldO1xuICAgICAgICAgICAgICAgICAgICB2YXIgZGFzR3JvdXAgPSBuZXcgREFTR3JvdXAoKTtcbiAgICAgICAgICAgICAgICAgICAgZGFzR3JvdXAudHlwZSA9IGdyb3VwWE1MLmdldEF0dHJpYnV0ZSgndHlwZScpO1xuICAgICAgICAgICAgICAgICAgICBkYXNHcm91cC5pZCA9IGdyb3VwWE1MLmdldEF0dHJpYnV0ZSgnaWQnKTtcbiAgICAgICAgICAgICAgICAgICAgZGFzR3JvdXAubGlua3MgPSBkYXNMaW5rc09mKGdyb3VwWE1MKTtcbiAgICAgICAgICAgICAgICAgICAgZGFzR3JvdXAubm90ZXMgPSBkYXNOb3Rlc09mKGdyb3VwWE1MKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKCFkYXNGZWF0dXJlLmdyb3Vwcykge1xuICAgICAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5ncm91cHMgPSBuZXcgQXJyYXkoZGFzR3JvdXApO1xuICAgICAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5ncm91cHMucHVzaChkYXNHcm91cCk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAvLyBNYWdpYyBub3Rlcy4gIENoZWNrIHdpdGggVEFEIGJlZm9yZSBjaGFuZ2luZyB0aGlzLlxuICAgICAgICAgICAgICAgIGlmIChkYXNGZWF0dXJlLm5vdGVzKSB7XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIG5pID0gMDsgbmkgPCBkYXNGZWF0dXJlLm5vdGVzLmxlbmd0aDsgKytuaSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIG4gPSBkYXNGZWF0dXJlLm5vdGVzW25pXTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChuLmluZGV4T2YoJ0dlbmVuYW1lPScpID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgZ2cgPSBuZXcgREFTR3JvdXAoKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBnZy50eXBlPSdnZW5lJztcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBnZy5pZCA9IG4uc3Vic3RyaW5nKDkpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmICghZGFzRmVhdHVyZS5ncm91cHMpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5ncm91cHMgPSBuZXcgQXJyYXkoZ2cpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUuZ3JvdXBzLnB1c2goZ2cpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBwZWMgPSBmZWF0dXJlLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdQQVJUJyk7XG4gICAgICAgICAgICAgICAgICAgIGlmIChwZWMubGVuZ3RoID4gMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHBhcnRzID0gW107XG4gICAgICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBwaSA9IDA7IHBpIDwgcGVjLmxlbmd0aDsgKytwaSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHBhcnRzLnB1c2gocGVjW3BpXS5nZXRBdHRyaWJ1dGUoJ2lkJykpO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5wYXJ0cyA9IHBhcnRzO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHBlYyA9IGZlYXR1cmUuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1BBUkVOVCcpO1xuICAgICAgICAgICAgICAgICAgICBpZiAocGVjLmxlbmd0aCA+IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBwYXJlbnRzID0gW107XG4gICAgICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBwaSA9IDA7IHBpIDwgcGVjLmxlbmd0aDsgKytwaSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHBhcmVudHMucHVzaChwZWNbcGldLmdldEF0dHJpYnV0ZSgnaWQnKSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLnBhcmVudHMgPSBwYXJlbnRzO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIGZlYXR1cmVzLnB1c2goZGFzRmVhdHVyZSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgICAgICAgICBcbiAgICAgICAgY2FsbGJhY2soZmVhdHVyZXMsIHVuZGVmaW5lZCwgc2VnbWVudE1hcCk7XG4gICAgfSxcbiAgICBmdW5jdGlvbiAoZXJyKSB7XG4gICAgICAgIGNhbGxiYWNrKFtdLCBlcnIpO1xuICAgIH0pO1xufVxuXG5mdW5jdGlvbiBEQVNBbGlnbm1lbnQodHlwZSkge1xuICAgIHRoaXMudHlwZSA9IHR5cGU7XG4gICAgdGhpcy5vYmplY3RzID0ge307XG4gICAgdGhpcy5ibG9ja3MgPSBbXTtcbn1cblxuREFTU291cmNlLnByb3RvdHlwZS5hbGlnbm1lbnRzID0gZnVuY3Rpb24oc2VnbWVudCwgb3B0aW9ucywgY2FsbGJhY2spIHtcbiAgICB2YXIgZGFzVVJJID0gdGhpcy5kYXNCYXNlVVJJICsgJ2FsaWdubWVudD9xdWVyeT0nICsgc2VnbWVudDtcbiAgICB0aGlzLmRvQ3Jvc3NEb21haW5SZXF1ZXN0KGRhc1VSSSwgZnVuY3Rpb24ocmVzcG9uc2VYTUwpIHtcbiAgICAgICAgaWYgKCFyZXNwb25zZVhNTCkge1xuICAgICAgICAgICAgY2FsbGJhY2soW10sICdGYWlsZWQgcmVxdWVzdCAnICsgZGFzVVJJKTtcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuXG4gICAgICAgIHZhciBhbGlnbm1lbnRzID0gW107XG4gICAgICAgIHZhciBhbGlYTUxzID0gcmVzcG9uc2VYTUwuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ2FsaWdubWVudCcpO1xuICAgICAgICBmb3IgKHZhciBhaSA9IDA7IGFpIDwgYWxpWE1Mcy5sZW5ndGg7ICsrYWkpIHtcbiAgICAgICAgICAgIHZhciBhbGlYTUwgPSBhbGlYTUxzW2FpXTtcbiAgICAgICAgICAgIHZhciBhbGkgPSBuZXcgREFTQWxpZ25tZW50KGFsaVhNTC5nZXRBdHRyaWJ1dGUoJ2FsaWduVHlwZScpKTtcbiAgICAgICAgICAgIHZhciBvYmpYTUxzID0gYWxpWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdhbGlnbk9iamVjdCcpO1xuICAgICAgICAgICAgZm9yICh2YXIgb2kgPSAwOyBvaSA8IG9ialhNTHMubGVuZ3RoOyArK29pKSB7XG4gICAgICAgICAgICAgICAgdmFyIG9ialhNTCA9IG9ialhNTHNbb2ldO1xuICAgICAgICAgICAgICAgIHZhciBvYmogPSB7XG4gICAgICAgICAgICAgICAgICAgIGlkOiAgICAgICAgICBvYmpYTUwuZ2V0QXR0cmlidXRlKCdpbnRPYmplY3RJZCcpLFxuICAgICAgICAgICAgICAgICAgICBhY2Nlc3Npb246ICAgb2JqWE1MLmdldEF0dHJpYnV0ZSgnZGJBY2Nlc3Npb25JZCcpLFxuICAgICAgICAgICAgICAgICAgICB2ZXJzaW9uOiAgICAgb2JqWE1MLmdldEF0dHJpYnV0ZSgnb2JqZWN0VmVyc2lvbicpLFxuICAgICAgICAgICAgICAgICAgICBkYlNvdXJjZTogICAgb2JqWE1MLmdldEF0dHJpYnV0ZSgnZGJTb3VyY2UnKSxcbiAgICAgICAgICAgICAgICAgICAgZGJWZXJzaW9uOiAgIG9ialhNTC5nZXRBdHRyaWJ1dGUoJ2RiVmVyc2lvbicpXG4gICAgICAgICAgICAgICAgfTtcbiAgICAgICAgICAgICAgICBhbGkub2JqZWN0c1tvYmouaWRdID0gb2JqO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgXG4gICAgICAgICAgICB2YXIgYmxvY2tYTUxzID0gYWxpWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdibG9jaycpO1xuICAgICAgICAgICAgZm9yICh2YXIgYmkgPSAwOyBiaSA8IGJsb2NrWE1Mcy5sZW5ndGg7ICsrYmkpIHtcbiAgICAgICAgICAgICAgICB2YXIgYmxvY2tYTUwgPSBibG9ja1hNTHNbYmldO1xuICAgICAgICAgICAgICAgIHZhciBibG9jayA9IHtcbiAgICAgICAgICAgICAgICAgICAgb3JkZXI6ICAgICAgYmxvY2tYTUwuZ2V0QXR0cmlidXRlKCdibG9ja09yZGVyJyksXG4gICAgICAgICAgICAgICAgICAgIHNlZ21lbnRzOiAgIFtdXG4gICAgICAgICAgICAgICAgfTtcbiAgICAgICAgICAgICAgICB2YXIgc2VnWE1McyA9IGJsb2NrWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdzZWdtZW50Jyk7XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgc2kgPSAwOyBzaSA8IHNlZ1hNTHMubGVuZ3RoOyArK3NpKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWdYTUwgPSBzZWdYTUxzW3NpXTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHNlZyA9IHtcbiAgICAgICAgICAgICAgICAgICAgICAgIG9iamVjdDogICAgICBzZWdYTUwuZ2V0QXR0cmlidXRlKCdpbnRPYmplY3RJZCcpLFxuICAgICAgICAgICAgICAgICAgICAgICAgbWluOiAgICAgICAgIHNlZ1hNTC5nZXRBdHRyaWJ1dGUoJ3N0YXJ0JyksXG4gICAgICAgICAgICAgICAgICAgICAgICBtYXg6ICAgICAgICAgc2VnWE1MLmdldEF0dHJpYnV0ZSgnZW5kJyksXG4gICAgICAgICAgICAgICAgICAgICAgICBzdHJhbmQ6ICAgICAgc2VnWE1MLmdldEF0dHJpYnV0ZSgnc3RyYW5kJyksXG4gICAgICAgICAgICAgICAgICAgICAgICBjaWdhcjogICAgICAgZWxlbWVudFZhbHVlKHNlZ1hNTCwgJ2NpZ2FyJylcbiAgICAgICAgICAgICAgICAgICAgfTtcbiAgICAgICAgICAgICAgICAgICAgYmxvY2suc2VnbWVudHMucHVzaChzZWcpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBhbGkuYmxvY2tzLnB1c2goYmxvY2spO1xuICAgICAgICAgICAgfSAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICBhbGlnbm1lbnRzLnB1c2goYWxpKTtcbiAgICAgICAgfVxuICAgICAgICBjYWxsYmFjayhhbGlnbm1lbnRzKTtcbiAgICB9KTtcbn1cblxuXG5mdW5jdGlvbiBEQVNTdHlsZXNoZWV0KCkge1xuLypcbiAgICB0aGlzLmhpZ2hab29tU3R5bGVzID0gbmV3IE9iamVjdCgpO1xuICAgIHRoaXMubWVkaXVtWm9vbVN0eWxlcyA9IG5ldyBPYmplY3QoKTtcbiAgICB0aGlzLmxvd1pvb21TdHlsZXMgPSBuZXcgT2JqZWN0KCk7XG4qL1xuXG4gICAgdGhpcy5zdHlsZXMgPSBbXTtcbn1cblxuREFTU3R5bGVzaGVldC5wcm90b3R5cGUucHVzaFN0eWxlID0gZnVuY3Rpb24oZmlsdGVycywgem9vbSwgc3R5bGUpIHtcbiAgICAvKlxuXG4gICAgaWYgKCF6b29tKSB7XG4gICAgICAgIHRoaXMuaGlnaFpvb21TdHlsZXNbdHlwZV0gPSBzdHlsZTtcbiAgICAgICAgdGhpcy5tZWRpdW1ab29tU3R5bGVzW3R5cGVdID0gc3R5bGU7XG4gICAgICAgIHRoaXMubG93Wm9vbVN0eWxlc1t0eXBlXSA9IHN0eWxlO1xuICAgIH0gZWxzZSBpZiAoem9vbSA9PSAnaGlnaCcpIHtcbiAgICAgICAgdGhpcy5oaWdoWm9vbVN0eWxlc1t0eXBlXSA9IHN0eWxlO1xuICAgIH0gZWxzZSBpZiAoem9vbSA9PSAnbWVkaXVtJykge1xuICAgICAgICB0aGlzLm1lZGl1bVpvb21TdHlsZXNbdHlwZV0gPSBzdHlsZTtcbiAgICB9IGVsc2UgaWYgKHpvb20gPT0gJ2xvdycpIHtcbiAgICAgICAgdGhpcy5sb3dab29tU3R5bGVzW3R5cGVdID0gc3R5bGU7XG4gICAgfVxuXG4gICAgKi9cblxuICAgIGlmICghZmlsdGVycykge1xuICAgICAgICBmaWx0ZXJzID0ge3R5cGU6ICdkZWZhdWx0J307XG4gICAgfVxuICAgIHZhciBzdHlsZUhvbGRlciA9IHNoYWxsb3dDb3B5KGZpbHRlcnMpO1xuICAgIGlmICh6b29tKSB7XG4gICAgICAgIHN0eWxlSG9sZGVyLnpvb20gPSB6b29tO1xuICAgIH1cbiAgICBzdHlsZUhvbGRlci5zdHlsZSA9IHN0eWxlO1xuICAgIHRoaXMuc3R5bGVzLnB1c2goc3R5bGVIb2xkZXIpO1xufVxuXG5mdW5jdGlvbiBEQVNTdHlsZSgpIHtcbn1cblxuZnVuY3Rpb24gcGFyc2VHcmFkaWVudChncmFkKSB7XG4gICAgdmFyIHN0ZXBzID0gZ3JhZC5nZXRBdHRyaWJ1dGUoJ3N0ZXBzJyk7XG4gICAgaWYgKHN0ZXBzKSB7XG4gICAgICAgIHN0ZXBzID0gc3RlcHN8MDtcbiAgICB9IGVsc2Uge1xuICAgICAgICBzdGVwcyA9IDUwO1xuICAgIH1cblxuXG4gICAgdmFyIHN0b3BzID0gW107XG4gICAgdmFyIGNvbG9ycyA9IFtdO1xuICAgIHZhciBzZSA9IGdyYWQuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1NUT1AnKTtcbiAgICBmb3IgKHZhciBzaSA9IDA7IHNpIDwgc2UubGVuZ3RoOyArK3NpKSB7XG4gICAgICAgIHZhciBzdG9wID0gc2Vbc2ldO1xuICAgICAgICBzdG9wcy5wdXNoKDEuMCAqIHN0b3AuZ2V0QXR0cmlidXRlKCdzY29yZScpKTtcbiAgICAgICAgY29sb3JzLnB1c2goc3RvcC5maXJzdENoaWxkLm5vZGVWYWx1ZSk7XG4gICAgfVxuXG4gICAgcmV0dXJuIG1ha2VDb2xvdXJTdGVwcyhzdGVwcywgc3RvcHMsIGNvbG9ycyk7XG59XG5cbkRBU1NvdXJjZS5wcm90b3R5cGUuc3R5bGVzaGVldCA9IGZ1bmN0aW9uKHN1Y2Nlc3NDQiwgZmFpbHVyZUNCKSB7XG4gICAgdmFyIGRhc1VSSSwgY3JlZHMgPSB0aGlzLmNyZWRlbnRpYWxzO1xuICAgIGlmICh0aGlzLnN0eWxlc2hlZXRfdXJpKSB7XG4gICAgICAgIGRhc1VSSSA9IHRoaXMuc3R5bGVzaGVldF91cmk7XG4gICAgICAgIGNyZWRzID0gZmFsc2U7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgZGFzVVJJID0gdGhpcy5kYXNCYXNlVVJJICsgJ3N0eWxlc2hlZXQnO1xuICAgIH1cblxuICAgIGRvQ3Jvc3NEb21haW5SZXF1ZXN0KGRhc1VSSSwgZnVuY3Rpb24ocmVzcG9uc2VYTUwpIHtcbiAgICAgICAgaWYgKCFyZXNwb25zZVhNTCkge1xuICAgICAgICAgICAgaWYgKGZhaWx1cmVDQikge1xuICAgICAgICAgICAgICAgIGZhaWx1cmVDQigpO1xuICAgICAgICAgICAgfSBcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuICAgICAgICB2YXIgc3R5bGVzaGVldCA9IG5ldyBEQVNTdHlsZXNoZWV0KCk7XG4gICAgICAgIHZhciB0eXBlWE1McyA9IHJlc3BvbnNlWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdUWVBFJyk7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdHlwZVhNTHMubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgIHZhciB0eXBlU3R5bGUgPSB0eXBlWE1Mc1tpXTtcbiAgICAgICAgICAgIFxuICAgICAgICAgICAgdmFyIGZpbHRlciA9IHt9O1xuICAgICAgICAgICAgZmlsdGVyLnR5cGUgPSB0eXBlU3R5bGUuZ2V0QXR0cmlidXRlKCdpZCcpOyAvLyBBbSBJIHJpZ2h0IGluIHRoaW5raW5nIHRoYXQgdGhpcyBtYWtlcyBEQVNTVFlMRSBYTUwgaW52YWxpZD8gIFVnaC5cbiAgICAgICAgICAgIGZpbHRlci5sYWJlbCA9IHR5cGVTdHlsZS5nZXRBdHRyaWJ1dGUoJ2xhYmVsJyk7XG4gICAgICAgICAgICBmaWx0ZXIubWV0aG9kID0gdHlwZVN0eWxlLmdldEF0dHJpYnV0ZSgnbWV0aG9kJyk7XG4gICAgICAgICAgICB2YXIgZ2x5cGhYTUxzID0gdHlwZVN0eWxlLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdHTFlQSCcpO1xuICAgICAgICAgICAgZm9yICh2YXIgZ2kgPSAwOyBnaSA8IGdseXBoWE1Mcy5sZW5ndGg7ICsrZ2kpIHtcbiAgICAgICAgICAgICAgICB2YXIgZ2x5cGhYTUwgPSBnbHlwaFhNTHNbZ2ldO1xuICAgICAgICAgICAgICAgIHZhciB6b29tID0gZ2x5cGhYTUwuZ2V0QXR0cmlidXRlKCd6b29tJyk7XG4gICAgICAgICAgICAgICAgdmFyIGdseXBoID0gY2hpbGRFbGVtZW50T2YoZ2x5cGhYTUwpO1xuICAgICAgICAgICAgICAgIHZhciBzdHlsZSA9IG5ldyBEQVNTdHlsZSgpO1xuICAgICAgICAgICAgICAgIHN0eWxlLmdseXBoID0gZ2x5cGgubG9jYWxOYW1lO1xuICAgICAgICAgICAgICAgIHZhciBjaGlsZCA9IGdseXBoLmZpcnN0Q2hpbGQ7XG4gICAgICAgIFxuICAgICAgICAgICAgICAgIHdoaWxlIChjaGlsZCkge1xuICAgICAgICAgICAgICAgICAgICBpZiAoY2hpbGQubm9kZVR5cGUgPT0gTm9kZS5FTEVNRU5UX05PREUpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIC8vIGFsZXJ0KGNoaWxkLmxvY2FsTmFtZSk7XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoY2hpbGQubG9jYWxOYW1lID09ICdCR0dSQUQnKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgc3R5bGVbY2hpbGQubG9jYWxOYW1lXSA9IHBhcnNlR3JhZGllbnQoY2hpbGQpO1xuICAgICAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHsgICAgICBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBzdHlsZVtjaGlsZC5sb2NhbE5hbWVdID0gY2hpbGQuZmlyc3RDaGlsZC5ub2RlVmFsdWU7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgY2hpbGQgPSBjaGlsZC5uZXh0U2libGluZztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgc3R5bGVzaGVldC5wdXNoU3R5bGUoZmlsdGVyLCB6b29tLCBzdHlsZSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgc3VjY2Vzc0NCKHN0eWxlc2hlZXQpO1xuICAgIH0sIGNyZWRzKTtcbn1cblxuLy9cbi8vIHNvdXJjZXMgY29tbWFuZFxuLy8gXG5cbmZ1bmN0aW9uIERBU1JlZ2lzdHJ5KHVyaSwgb3B0cylcbntcbiAgICBvcHRzID0gb3B0cyB8fCB7fTtcbiAgICB0aGlzLnVyaSA9IHVyaTtcbiAgICB0aGlzLm9wdHMgPSBvcHRzOyAgIFxufVxuXG5EQVNSZWdpc3RyeS5wcm90b3R5cGUuc291cmNlcyA9IGZ1bmN0aW9uKGNhbGxiYWNrLCBmYWlsdXJlLCBvcHRzKVxue1xuICAgIGlmICghb3B0cykge1xuICAgICAgICBvcHRzID0ge307XG4gICAgfVxuXG4gICAgdmFyIGZpbHRlcnMgPSBbXTtcbiAgICBpZiAob3B0cy50YXhvbikge1xuICAgICAgICBmaWx0ZXJzLnB1c2goJ29yZ2FuaXNtPScgKyBvcHRzLnRheG9uKTtcbiAgICB9XG4gICAgaWYgKG9wdHMuYXV0aCkge1xuICAgICAgICBmaWx0ZXJzLnB1c2goJ2F1dGhvcml0eT0nICsgb3B0cy5hdXRoKTtcbiAgICB9XG4gICAgaWYgKG9wdHMudmVyc2lvbikge1xuICAgICAgICBmaWx0ZXJzLnB1c2goJ3ZlcnNpb249JyArIG9wdHMudmVyc2lvbik7XG4gICAgfVxuICAgIHZhciBxdXJpID0gdGhpcy51cmk7XG4gICAgaWYgKGZpbHRlcnMubGVuZ3RoID4gMCkge1xuICAgICAgICBxdXJpID0gcXVyaSArICc/JyArIGZpbHRlcnMuam9pbignJicpOyAgIC8vICcmJyBhcyBhIHNlcGFyYXRvciB0byBoYWNrIGFyb3VuZCBkYXNyZWdpc3RyeS5vcmcgYnVnLlxuICAgIH1cblxuICAgIGRvQ3Jvc3NEb21haW5SZXF1ZXN0KHF1cmksIGZ1bmN0aW9uKHJlc3BvbnNlWE1MKSB7XG4gICAgICAgIGlmICghcmVzcG9uc2VYTUwgJiYgZmFpbHVyZSkge1xuICAgICAgICAgICAgZmFpbHVyZSgpO1xuICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICB9XG5cbiAgICAgICAgdmFyIHNvdXJjZXMgPSBbXTsgICAgICAgXG4gICAgICAgIHZhciBzb3VyY2VYTUxzID0gcmVzcG9uc2VYTUwuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1NPVVJDRScpO1xuICAgICAgICBmb3IgKHZhciBzaSA9IDA7IHNpIDwgc291cmNlWE1Mcy5sZW5ndGg7ICsrc2kpIHtcbiAgICAgICAgICAgIHZhciBzb3VyY2VYTUwgPSBzb3VyY2VYTUxzW3NpXTtcbiAgICAgICAgICAgIHZhciB2ZXJzaW9uWE1McyA9IHNvdXJjZVhNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnVkVSU0lPTicpO1xuICAgICAgICAgICAgaWYgKHZlcnNpb25YTUxzLmxlbmd0aCA8IDEpIHtcbiAgICAgICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHZhciB2ZXJzaW9uWE1MID0gdmVyc2lvblhNTHNbMF07XG5cbiAgICAgICAgICAgIHZhciBjb29yZFhNTHMgPSB2ZXJzaW9uWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdDT09SRElOQVRFUycpO1xuICAgICAgICAgICAgdmFyIGNvb3JkcyA9IFtdO1xuICAgICAgICAgICAgZm9yICh2YXIgY2kgPSAwOyBjaSA8IGNvb3JkWE1Mcy5sZW5ndGg7ICsrY2kpIHtcbiAgICAgICAgICAgICAgICB2YXIgY29vcmRYTUwgPSBjb29yZFhNTHNbY2ldO1xuICAgICAgICAgICAgICAgIHZhciBjb29yZCA9IG5ldyBEQVNDb29yZHMoKTtcbiAgICAgICAgICAgICAgICBjb29yZC5hdXRoID0gY29vcmRYTUwuZ2V0QXR0cmlidXRlKCdhdXRob3JpdHknKTtcbiAgICAgICAgICAgICAgICBjb29yZC50YXhvbiA9IGNvb3JkWE1MLmdldEF0dHJpYnV0ZSgndGF4aWQnKTtcbiAgICAgICAgICAgICAgICBjb29yZC52ZXJzaW9uID0gY29vcmRYTUwuZ2V0QXR0cmlidXRlKCd2ZXJzaW9uJyk7XG4gICAgICAgICAgICAgICAgY29vcmRzLnB1c2goY29vcmQpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgXG4gICAgICAgICAgICB2YXIgY2FwcyA9IFtdO1xuICAgICAgICAgICAgdmFyIGNhcFhNTHMgPSB2ZXJzaW9uWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdDQVBBQklMSVRZJyk7XG4gICAgICAgICAgICB2YXIgdXJpO1xuICAgICAgICAgICAgZm9yICh2YXIgY2kgPSAwOyBjaSA8IGNhcFhNTHMubGVuZ3RoOyArK2NpKSB7XG4gICAgICAgICAgICAgICAgdmFyIGNhcFhNTCA9IGNhcFhNTHNbY2ldO1xuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIGNhcHMucHVzaChjYXBYTUwuZ2V0QXR0cmlidXRlKCd0eXBlJykpO1xuXG4gICAgICAgICAgICAgICAgaWYgKGNhcFhNTC5nZXRBdHRyaWJ1dGUoJ3R5cGUnKSA9PSAnZGFzMTpmZWF0dXJlcycpIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGZlcCA9IGNhcFhNTC5nZXRBdHRyaWJ1dGUoJ3F1ZXJ5X3VyaScpO1xuICAgICAgICAgICAgICAgICAgICB1cmkgPSBmZXAuc3Vic3RyaW5nKDAsIGZlcC5sZW5ndGggLSAoJ2ZlYXR1cmVzJy5sZW5ndGgpKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHZhciBwcm9wcyA9IHt9O1xuICAgICAgICAgICAgdmFyIHByb3BYTUxzID0gdmVyc2lvblhNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnUFJPUCcpO1xuICAgICAgICAgICAgZm9yICh2YXIgcGkgPSAwOyBwaSA8IHByb3BYTUxzLmxlbmd0aDsgKytwaSkge1xuICAgICAgICAgICAgICAgIHB1c2hvKHByb3BzLCBwcm9wWE1Mc1twaV0uZ2V0QXR0cmlidXRlKCduYW1lJyksIHByb3BYTUxzW3BpXS5nZXRBdHRyaWJ1dGUoJ3ZhbHVlJykpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgXG4gICAgICAgICAgICBpZiAodXJpKSB7XG4gICAgICAgICAgICAgICAgdmFyIHNvdXJjZSA9IG5ldyBEQVNTb3VyY2UodXJpLCB7XG4gICAgICAgICAgICAgICAgICAgIHNvdXJjZV91cmk6IHNvdXJjZVhNTC5nZXRBdHRyaWJ1dGUoJ3VyaScpLFxuICAgICAgICAgICAgICAgICAgICBuYW1lOiAgc291cmNlWE1MLmdldEF0dHJpYnV0ZSgndGl0bGUnKSxcbiAgICAgICAgICAgICAgICAgICAgZGVzYzogIHNvdXJjZVhNTC5nZXRBdHRyaWJ1dGUoJ2Rlc2NyaXB0aW9uJyksXG4gICAgICAgICAgICAgICAgICAgIGNvb3JkczogY29vcmRzLFxuICAgICAgICAgICAgICAgICAgICBwcm9wczogcHJvcHMsXG4gICAgICAgICAgICAgICAgICAgIGNhcGFiaWxpdGllczogY2Fwc1xuICAgICAgICAgICAgICAgIH0pO1xuICAgICAgICAgICAgICAgIHNvdXJjZXMucHVzaChzb3VyY2UpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIFxuICAgICAgICBjYWxsYmFjayhzb3VyY2VzKTtcbiAgICB9KTtcbn1cblxuXG4vL1xuLy8gVXRpbGl0eSBmdW5jdGlvbnNcbi8vXG5cbmZ1bmN0aW9uIGVsZW1lbnRWYWx1ZShlbGVtZW50LCB0YWcpXG57XG4gICAgdmFyIGNoaWxkcmVuID0gZWxlbWVudC5nZXRFbGVtZW50c0J5VGFnTmFtZSh0YWcpO1xuICAgIGlmIChjaGlsZHJlbi5sZW5ndGggPiAwICYmIGNoaWxkcmVuWzBdLmZpcnN0Q2hpbGQpIHtcbiAgICAgICAgdmFyIGMgPSBjaGlsZHJlblswXTtcbiAgICAgICAgaWYgKGMuY2hpbGROb2Rlcy5sZW5ndGggPT0gMSkge1xuICAgICAgICAgICAgcmV0dXJuIGMuZmlyc3RDaGlsZC5ub2RlVmFsdWU7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB2YXIgcyA9ICcnO1xuICAgICAgICAgICAgZm9yICh2YXIgbmkgPSAwOyBuaSA8IGMuY2hpbGROb2Rlcy5sZW5ndGg7ICsrbmkpIHtcbiAgICAgICAgICAgICAgICBzICs9IGMuY2hpbGROb2Rlc1tuaV0ubm9kZVZhbHVlO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuIHM7XG4gICAgICAgIH1cblxuICAgIH0gZWxzZSB7XG4gICAgICAgIHJldHVybiBudWxsO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gY2hpbGRFbGVtZW50T2YoZWxlbWVudClcbntcbiAgICBpZiAoZWxlbWVudC5oYXNDaGlsZE5vZGVzKCkpIHtcbiAgICAgICAgdmFyIGNoaWxkID0gZWxlbWVudC5maXJzdENoaWxkO1xuICAgICAgICBkbyB7XG4gICAgICAgICAgICBpZiAoY2hpbGQubm9kZVR5cGUgPT0gTm9kZS5FTEVNRU5UX05PREUpIHtcbiAgICAgICAgICAgICAgICByZXR1cm4gY2hpbGQ7XG4gICAgICAgICAgICB9IFxuICAgICAgICAgICAgY2hpbGQgPSBjaGlsZC5uZXh0U2libGluZztcbiAgICAgICAgfSB3aGlsZSAoY2hpbGQgIT0gbnVsbCk7XG4gICAgfVxuICAgIHJldHVybiBudWxsO1xufVxuXG5cbmZ1bmN0aW9uIGRhc0xpbmtzT2YoZWxlbWVudClcbntcbiAgICB2YXIgbGlua3MgPSBuZXcgQXJyYXkoKTtcbiAgICB2YXIgbWF5YmVMaW5rQ2hpbGRlbiA9IGVsZW1lbnQuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ0xJTksnKTtcbiAgICBmb3IgKHZhciBjaSA9IDA7IGNpIDwgbWF5YmVMaW5rQ2hpbGRlbi5sZW5ndGg7ICsrY2kpIHtcbiAgICAgICAgdmFyIGxpbmtYTUwgPSBtYXliZUxpbmtDaGlsZGVuW2NpXTtcbiAgICAgICAgaWYgKGxpbmtYTUwucGFyZW50Tm9kZSA9PSBlbGVtZW50KSB7XG4gICAgICAgICAgICBsaW5rcy5wdXNoKG5ldyBEQVNMaW5rKGxpbmtYTUwuZmlyc3RDaGlsZCA/IGxpbmtYTUwuZmlyc3RDaGlsZC5ub2RlVmFsdWUgOiAnVW5rbm93bicsIGxpbmtYTUwuZ2V0QXR0cmlidXRlKCdocmVmJykpKTtcbiAgICAgICAgfVxuICAgIH1cbiAgICBcbiAgICByZXR1cm4gbGlua3M7XG59XG5cbmZ1bmN0aW9uIGRhc05vdGVzT2YoZWxlbWVudClcbntcbiAgICB2YXIgbm90ZXMgPSBbXTtcbiAgICB2YXIgbWF5YmVOb3RlcyA9IGVsZW1lbnQuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ05PVEUnKTtcbiAgICBmb3IgKHZhciBuaSA9IDA7IG5pIDwgbWF5YmVOb3Rlcy5sZW5ndGg7ICsrbmkpIHtcbiAgICAgICAgaWYgKG1heWJlTm90ZXNbbmldLmZpcnN0Q2hpbGQpIHtcbiAgICAgICAgICAgIG5vdGVzLnB1c2gobWF5YmVOb3Rlc1tuaV0uZmlyc3RDaGlsZC5ub2RlVmFsdWUpO1xuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiBub3Rlcztcbn1cblxuZnVuY3Rpb24gZG9Dcm9zc0RvbWFpblJlcXVlc3QodXJsLCBoYW5kbGVyLCBjcmVkZW50aWFscywgY3VzdEF1dGgpIHtcbiAgICAvLyBUT0RPOiBleHBsaWNpdCBlcnJvciBoYW5kbGVycz9cblxuICAgIGlmICh3aW5kb3cuWERvbWFpblJlcXVlc3QpIHtcbiAgICAgICAgdmFyIHJlcSA9IG5ldyBYRG9tYWluUmVxdWVzdCgpO1xuICAgICAgICByZXEub25sb2FkID0gZnVuY3Rpb24oKSB7XG4gICAgICAgICAgICB2YXIgZG9tID0gbmV3IEFjdGl2ZVhPYmplY3QoXCJNaWNyb3NvZnQuWE1MRE9NXCIpO1xuICAgICAgICAgICAgZG9tLmFzeW5jID0gZmFsc2U7XG4gICAgICAgICAgICBkb20ubG9hZFhNTChyZXEucmVzcG9uc2VUZXh0KTtcbiAgICAgICAgICAgIGhhbmRsZXIoZG9tKTtcbiAgICAgICAgfVxuICAgICAgICByZXEub3BlbihcImdldFwiLCB1cmwpO1xuICAgICAgICByZXEuc2VuZCgnJyk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdmFyIHJlcVN0YXJ0ID0gRGF0ZS5ub3coKTtcbiAgICAgICAgdmFyIHJlcSA9IG5ldyBYTUxIdHRwUmVxdWVzdCgpO1xuXG4gICAgICAgIHJlcS5vbnJlYWR5c3RhdGVjaGFuZ2UgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgICAgIGlmIChyZXEucmVhZHlTdGF0ZSA9PSA0KSB7XG4gICAgICAgICAgICAgIGlmIChyZXEuc3RhdHVzID49IDIwMCB8fCByZXEuc3RhdHVzID09IDApIHtcbiAgICAgICAgICAgICAgICAgIGhhbmRsZXIocmVxLnJlc3BvbnNlWE1MLCByZXEpO1xuICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH07XG4gICAgICAgIHJlcS5vcGVuKFwiZ2V0XCIsIHVybCwgdHJ1ZSk7XG4gICAgICAgIGlmIChjcmVkZW50aWFscykge1xuICAgICAgICAgICAgcmVxLndpdGhDcmVkZW50aWFscyA9IHRydWU7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKGN1c3RBdXRoKSB7XG4gICAgICAgICAgICByZXEuc2V0UmVxdWVzdEhlYWRlcignWC1EQVMtQXV0aG9yaXNhdGlvbicsIGN1c3RBdXRoKTtcbiAgICAgICAgfVxuICAgICAgICByZXEub3ZlcnJpZGVNaW1lVHlwZSgndGV4dC94bWwnKTtcbiAgICAgICAgcmVxLnNldFJlcXVlc3RIZWFkZXIoJ0FjY2VwdCcsICdhcHBsaWNhdGlvbi94bWwsKi8qJyk7XG4gICAgICAgIHJlcS5zZW5kKCcnKTtcbiAgICB9XG59XG5cbkRBU1NvdXJjZS5wcm90b3R5cGUuZG9Dcm9zc0RvbWFpblJlcXVlc3QgPSBmdW5jdGlvbih1cmwsIGhhbmRsZXIsIGVyckhhbmRsZXIpIHtcbiAgICB2YXIgY3VzdEF1dGg7XG4gICAgaWYgKHRoaXMueFVzZXIpIHtcbiAgICAgICAgY3VzdEF1dGggPSAnQmFzaWMgJyArIGJ0b2EodGhpcy54VXNlciArICc6JyArIHRoaXMueFBhc3MpO1xuICAgIH1cblxuICAgIHRyeSB7XG4gICAgICAgIHJldHVybiBkb0Nyb3NzRG9tYWluUmVxdWVzdCh1cmwsIGhhbmRsZXIsIHRoaXMuY3JlZGVudGlhbHMsIGN1c3RBdXRoKTtcbiAgICB9IGNhdGNoIChlcnIpIHtcbiAgICAgICAgaWYgKGVyckhhbmRsZXIpIHtcbiAgICAgICAgICAgIGVyckhhbmRsZXIoZXJyKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHRocm93IGVycjtcbiAgICAgICAgfVxuICAgIH1cbn1cblxuZnVuY3Rpb24gaXNEYXNCb29sZWFuVHJ1ZShzKSB7XG4gICAgcyA9ICgnJyArIHMpLnRvTG93ZXJDYXNlKCk7XG4gICAgcmV0dXJuIHM9PT0neWVzJyB8fCBzPT09J3RydWUnO1xufVxuXG5mdW5jdGlvbiBpc0Rhc0Jvb2xlYW5Ob3RGYWxzZShzKSB7XG4gICAgaWYgKCFzKVxuICAgICAgICByZXR1cm4gZmFsc2U7XG5cbiAgICBzID0gKCcnICsgcykudG9Mb3dlckNhc2UoKTtcbiAgICByZXR1cm4gcyE9PSdubycgfHwgcyE9PSdmYWxzZSc7XG59XG5cbmZ1bmN0aW9uIGNvcHlTdHlsZXNoZWV0KHNzKSB7XG4gICAgdmFyIG5zcyA9IHNoYWxsb3dDb3B5KHNzKTtcbiAgICBuc3Muc3R5bGVzID0gW107XG4gICAgZm9yICh2YXIgc2kgPSAwOyBzaSA8IHNzLnN0eWxlcy5sZW5ndGg7ICsrc2kpIHtcbiAgICAgICAgdmFyIHNoID0gbnNzLnN0eWxlc1tzaV0gPSBzaGFsbG93Q29weShzcy5zdHlsZXNbc2ldKTtcbiAgICAgICAgc2guX21ldGhvZFJFID0gc2guX2xhYmVsUkUgPSBzaC5fdHlwZVJFID0gdW5kZWZpbmVkO1xuICAgICAgICBzaC5zdHlsZSA9IHNoYWxsb3dDb3B5KHNoLnN0eWxlKTtcbiAgICAgICAgc2guc3R5bGUuaWQgPSB1bmRlZmluZWQ7XG4gICAgICAgIHNoLnN0eWxlLl9ncmFkaWVudCA9IHVuZGVmaW5lZDtcbiAgICB9XG4gICAgcmV0dXJuIG5zcztcbn1cblxuaWYgKHR5cGVvZihtb2R1bGUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIG1vZHVsZS5leHBvcnRzID0ge1xuICAgICAgICBEQVNHcm91cDogREFTR3JvdXAsXG4gICAgICAgIERBU0ZlYXR1cmU6IERBU0ZlYXR1cmUsXG4gICAgICAgIERBU1N0eWxlc2hlZXQ6IERBU1N0eWxlc2hlZXQsXG4gICAgICAgIERBU1N0eWxlOiBEQVNTdHlsZSxcbiAgICAgICAgREFTU291cmNlOiBEQVNTb3VyY2UsXG4gICAgICAgIERBU1NlZ21lbnQ6IERBU1NlZ21lbnQsXG4gICAgICAgIERBU1JlZ2lzdHJ5OiBEQVNSZWdpc3RyeSxcbiAgICAgICAgREFTU2VxdWVuY2U6IERBU1NlcXVlbmNlLFxuICAgICAgICBEQVNMaW5rOiBEQVNMaW5rLFxuXG4gICAgICAgIGlzRGFzQm9vbGVhblRydWU6IGlzRGFzQm9vbGVhblRydWUsXG4gICAgICAgIGlzRGFzQm9vbGVhbk5vdEZhbHNlOiBpc0Rhc0Jvb2xlYW5Ob3RGYWxzZSxcbiAgICAgICAgY29weVN0eWxlc2hlZXQ6IGNvcHlTdHlsZXNoZWV0XG4gICAgfTtcbn0iLCIoZnVuY3Rpb24gKGdsb2JhbCl7XG4vKiAtKi0gbW9kZTogamF2YXNjcmlwdDsgYy1iYXNpYy1vZmZzZXQ6IDQ7IGluZGVudC10YWJzLW1vZGU6IG5pbCAtKi0gKi9cblxuLy8gXG4vLyBEYWxsaWFuY2UgR2Vub21lIEV4cGxvcmVyXG4vLyAoYykgVGhvbWFzIERvd24gMjAwNi0yMDE0XG4vL1xuLy8gZmV0Y2h3b3JrZXIuanNcbi8vXG5cblwidXNlIHN0cmljdFwiO1xuXG52YXIgYmluID0gcmVxdWlyZSgnLi9iaW4nKTtcbnZhciBiYW0gPSByZXF1aXJlKCcuL2JhbScpO1xudmFyIGJpZ3dpZyA9IHJlcXVpcmUoJy4vYmlnd2lnJyk7XG5cbnZhciBjb25uZWN0aW9ucyA9IHt9O1xuXG52YXIgaWRTZWVkID0gMDtcblxuZ2xvYmFsLm5ld0lEID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuICdjbicgKyAoKytpZFNlZWQpO1xufVxuXG5wb3N0TWVzc2FnZSh7dGFnOiAnaW5pdCd9KTtcblxuc2VsZi5vbm1lc3NhZ2UgPSBmdW5jdGlvbihldmVudCkge1xuICAgIHZhciBkID0gZXZlbnQuZGF0YTtcbiAgICB2YXIgY29tbWFuZCA9IGV2ZW50LmRhdGEuY29tbWFuZDtcbiAgICB2YXIgdGFnID0gZXZlbnQuZGF0YS50YWc7XG5cbiAgICBpZiAoY29tbWFuZCA9PT0gJ2Nvbm5lY3RCQU0nKSB7XG4gICAgICAgIHZhciBpZCA9IG5ld0lEKCk7XG5cbiAgICAgICAgdmFyIGJhbUYsIGJhaUY7XG4gICAgICAgIGlmIChkLmJsb2IpIHtcbiAgICAgICAgICAgIGJhbUYgPSBuZXcgYmluLkJsb2JGZXRjaGFibGUoZC5ibG9iKTtcbiAgICAgICAgICAgIGJhaUYgPSBuZXcgYmluLkJsb2JGZXRjaGFibGUoZC5pbmRleEJsb2IpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgYmFtRiA9IG5ldyBiaW4uVVJMRmV0Y2hhYmxlKGQudXJpLCB7Y3JlZGVudGlhbHM6IGQuY3JlZGVudGlhbHN9KTtcbiAgICAgICAgICAgIGJhaUYgPSBuZXcgYmluLlVSTEZldGNoYWJsZShkLmluZGV4VXJpLCB7Y3JlZGVudGlhbHM6IGQuY3JlZGVudGlhbHN9KTtcbiAgICAgICAgfVxuXG4gICAgICAgIGJhbS5tYWtlQmFtKGJhbUYsIGJhaUYsIGZ1bmN0aW9uKGJhbU9iaiwgZXJyKSB7XG4gICAgICAgICAgICBpZiAoYmFtT2JqKSB7XG4gICAgICAgICAgICAgICAgY29ubmVjdGlvbnNbaWRdID0gbmV3IEJBTVdvcmtlckZldGNoZXIoYmFtT2JqKTtcbiAgICAgICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogaWR9KTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCBlcnJvcjogZXJyIHx8IFwiQ291bGRuJ3QgZmV0Y2ggQkFNXCJ9KTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSk7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAnY29ubmVjdEJCSScpIHtcbiAgICAgICAgdmFyIGlkID0gbmV3SUQoKTtcbiAgICAgICAgdmFyIGJiaTtcbiAgICAgICAgaWYgKGQuYmxvYikge1xuICAgICAgICAgICAgYmJpID0gbmV3IGJpbi5CbG9iRmV0Y2hhYmxlKGQuYmxvYik7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBiYmkgPSBuZXcgYmluLlVSTEZldGNoYWJsZShkLnVyaSwge2NyZWRlbnRpYWxzOiBkLmNyZWRlbnRpYWxzfSk7XG4gICAgICAgIH1cblxuICAgICAgICBiaWd3aWcubWFrZUJ3ZyhiYmksIGZ1bmN0aW9uKGJ3ZywgZXJyKSB7XG4gICAgICAgICAgICBpZiAoYndnKSB7XG4gICAgICAgICAgICAgICAgY29ubmVjdGlvbnNbaWRdID0gbmV3IEJCSVdvcmtlckZldGNoZXIoYndnKTtcbiAgICAgICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogaWR9KTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCBlcnJvcjogZXJyIHx8IFwiQ291bGRuJ3QgZmV0Y2ggQkJJXCJ9KTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSwgZC51cmkpO1xuICAgIH0gZWxzZSBpZiAoY29tbWFuZCA9PT0gJ2ZldGNoJykge1xuICAgICAgICB2YXIgY29uID0gY29ubmVjdGlvbnNbZXZlbnQuZGF0YS5jb25uZWN0aW9uXTtcbiAgICAgICAgaWYgKCFjb24pIHtcbiAgICAgICAgICAgIHJldHVybiBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiAnTm8gc3VjaCBjb25uZWN0aW9uOiAnICsgZXZlbnQuZGF0YS5jb25uZWN0aW9ufSk7XG4gICAgICAgIH1cblxuICAgICAgICBjb24uZmV0Y2goZC50YWcsIGQuY2hyLCBkLm1pbiwgZC5tYXgsIGQuem9vbSwgZC5vcHRzKTtcbiAgICB9IGVsc2UgaWYgKGNvbW1hbmQgPT09ICdsZWFwJykge1xuICAgICAgICB2YXIgY29uID0gY29ubmVjdGlvbnNbZXZlbnQuZGF0YS5jb25uZWN0aW9uXTtcbiAgICAgICAgaWYgKCFjb24pIHtcbiAgICAgICAgICAgIHJldHVybiBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiAnTm8gc3VjaCBjb25uZWN0aW9uOiAnICsgZXZlbnQuZGF0YS5jb25uZWN0aW9ufSk7XG4gICAgICAgIH1cblxuICAgICAgICBjb24ubGVhcChkLnRhZywgZC5jaHIsIGQucG9zLCBkLmRpcik7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAncXVhbnRMZWFwJykge1xuICAgICAgICB2YXIgY29uID0gY29ubmVjdGlvbnNbZXZlbnQuZGF0YS5jb25uZWN0aW9uXTtcbiAgICAgICAgaWYgKCFjb24pIHtcbiAgICAgICAgICAgIHJldHVybiBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiAnTm8gc3VjaCBjb25uZWN0aW9uOiAnICsgZXZlbnQuZGF0YS5jb25uZWN0aW9ufSk7XG4gICAgICAgIH1cblxuICAgICAgICBjb24ucXVhbnRMZWFwKGQudGFnLCBkLmNociwgZC5wb3MsIGQuZGlyLCBkLnRocmVzaG9sZCwgZC51bmRlcik7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAnbWV0YScpIHtcbiAgICAgICAgdmFyIGNvbiA9IGNvbm5lY3Rpb25zW2V2ZW50LmRhdGEuY29ubmVjdGlvbl07XG4gICAgICAgIGlmICghY29uKSB7XG4gICAgICAgICAgICByZXR1cm4gcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCBlcnJvcjogJ05vIHN1Y2ggY29ubmVjdGlvbjogJyArIGV2ZW50LmRhdGEuY29ubmVjdGlvbn0pO1xuICAgICAgICB9XG5cbiAgICAgICAgY29uLm1ldGEoZC50YWcpO1xuICAgIH0gZWxzZSBpZiAoY29tbWFuZCA9PT0gJ3NlYXJjaCcpIHtcbiAgICAgICAgdmFyIGNvbiA9IGNvbm5lY3Rpb25zW2V2ZW50LmRhdGEuY29ubmVjdGlvbl07XG4gICAgICAgIGlmICghY29uKSB7XG4gICAgICAgICAgICByZXR1cm4gcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCBlcnJvcjogJ05vIHN1Y2ggY29ubmVjdGlvbjogJyArIGV2ZW50LmRhdGEuY29ubmVjdGlvbn0pO1xuICAgICAgICB9XG5cbiAgICAgICAgY29uLnNlYXJjaChkLnRhZywgZC5xdWVyeSwgZC5pbmRleCk7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAnZGF0ZScpIHtcbiAgICAgICAgcmV0dXJuIHBvc3RNZXNzYWdlKHt0YWc6IHRhZywgcmVzdWx0OiBEYXRlLm5vdygpfDB9KTtcbiAgICB9IGVsc2Uge1xuICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiAnQmFkIGNvbW1hbmQgJyArIGNvbW1hbmR9KTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIEJBTVdvcmtlckZldGNoZXIoYmFtKSB7XG4gICAgdGhpcy5iYW0gPSBiYW07XG59XG5cbkJBTVdvcmtlckZldGNoZXIucHJvdG90eXBlLmZldGNoID0gZnVuY3Rpb24odGFnLCBjaHIsIG1pbiwgbWF4LCB6b29tLCBvcHRzKSB7XG4gICAgb3B0cyA9IG9wdHMgfHwge307XG4gICAgdGhpcy5iYW0uZmV0Y2goY2hyLCBtaW4sIG1heCwgZnVuY3Rpb24ocmVjb3JkcywgZXJyKSB7XG4gICAgICAgIGlmIChyZWNvcmRzKSB7XG4gICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogcmVjb3JkcywgdGltZTogRGF0ZS5ub3coKXwwfSk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiBlcnJ9KTtcbiAgICAgICAgfVxuICAgIH0sIG9wdHMpO1xufVxuXG5mdW5jdGlvbiBCQklXb3JrZXJGZXRjaGVyKGJiaSkge1xuICAgIHRoaXMuYmJpID0gYmJpO1xufVxuXG5CQklXb3JrZXJGZXRjaGVyLnByb3RvdHlwZS5mZXRjaCA9IGZ1bmN0aW9uKHRhZywgY2hyLCBtaW4sIG1heCwgem9vbSkge1xuICAgIGlmICh0eXBlb2Yoem9vbSkgIT09ICdudW1iZXInKVxuICAgICAgICB6b29tID0gLTE7XG5cbiAgICB2YXIgZGF0YTtcbiAgICBpZiAoem9vbSA8IDApIHtcbiAgICAgICAgZGF0YSA9IHRoaXMuYmJpLmdldFVuem9vbWVkVmlldygpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIGRhdGEgPSB0aGlzLmJiaS5nZXRab29tZWRWaWV3KHpvb20pO1xuICAgIH1cblxuICAgIGRhdGEucmVhZFdpZ0RhdGEoY2hyLCBtaW4sIG1heCwgZnVuY3Rpb24oZmVhdHVyZXMpIHtcbiAgICAgICAgcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCByZXN1bHQ6IGZlYXR1cmVzfSk7XG4gICAgfSk7XG59XG5cbkJCSVdvcmtlckZldGNoZXIucHJvdG90eXBlLm1ldGEgPSBmdW5jdGlvbih0YWcpIHtcbiAgICB2YXIgc2NhbGVzID0gWzFdO1xuICAgIGZvciAodmFyIHogPSAwOyB6IDwgdGhpcy5iYmkuem9vbUxldmVscy5sZW5ndGg7ICsreikge1xuICAgICAgICBzY2FsZXMucHVzaCh0aGlzLmJiaS56b29tTGV2ZWxzW3pdLnJlZHVjdGlvbik7XG4gICAgfVxuXG4gICAgdmFyIHRoaXNCID0gdGhpcztcbiAgICB2YXIgbWV0YSA9IHt0eXBlOiB0aGlzLmJiaS50eXBlLFxuICAgICAgICAgICAgICAgIHpvb21MZXZlbHM6IHNjYWxlcyxcbiAgICAgICAgICAgICAgICBmaWVsZENvdW50OiB0aGlzLmJiaS5maWVsZENvdW50LFxuICAgICAgICAgICAgICAgIGRlZmluZWRGaWVsZENvdW50OiB0aGlzLmJiaS5kZWZpbmVkRmllbGRDb3VudCxcbiAgICAgICAgICAgICAgICBzY2hlbWE6IHRoaXMuYmJpLnNjaGVtYX07XG4gICAgaWYgKHRoaXMuYmJpLnR5cGUgPT09ICdiaWdiZWQnKSB7XG4gICAgICAgIHRoaXMuYmJpLmdldEV4dHJhSW5kaWNlcyhmdW5jdGlvbihlaSkge1xuICAgICAgICAgICAgaWYgKGVpKSB7XG4gICAgICAgICAgICAgICAgdGhpc0IuZXh0cmFJbmRpY2VzID0gZWk7XG4gICAgICAgICAgICAgICAgbWV0YS5leHRyYUluZGljZXMgPSBlaS5tYXAoZnVuY3Rpb24oaSkge3JldHVybiBpLmZpZWxkfSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogbWV0YX0pO1xuICAgICAgICB9KTtcbiAgICB9IGVsc2Uge1xuICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogbWV0YX0pO1xuICAgIH1cbn1cblxuQkJJV29ya2VyRmV0Y2hlci5wcm90b3R5cGUubGVhcCA9IGZ1bmN0aW9uKHRhZywgY2hyLCBwb3MsIGRpcikge1xuICAgIHRoaXMuYmJpLmdldFVuem9vbWVkVmlldygpLmdldEZpcnN0QWRqYWNlbnQoY2hyLCBwb3MsIGRpciwgZnVuY3Rpb24ocmVzdWx0LCBlcnIpIHtcbiAgICAgICAgcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCByZXN1bHQ6IHJlc3VsdCwgZXJyb3I6IGVycn0pO1xuICAgIH0pO1xufVxuXG5CQklXb3JrZXJGZXRjaGVyLnByb3RvdHlwZS5xdWFudExlYXAgPSBmdW5jdGlvbih0YWcsIGNociwgcG9zLCBkaXIsIHRocmVzaG9sZCwgdW5kZXIpIHtcbiAgICB0aGlzLmJiaS50aHJlc2hvbGRTZWFyY2goY2hyLCBwb3MsIGRpciwgdGhyZXNob2xkLCBmdW5jdGlvbihyZXN1bHQsIGVycikge1xuICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogcmVzdWx0LCBlcnJvcjogZXJyfSk7XG4gICAgfSk7XG59XG5cbkJCSVdvcmtlckZldGNoZXIucHJvdG90eXBlLnNlYXJjaCA9IGZ1bmN0aW9uKHRhZywgcXVlcnksIGluZGV4KSB7XG4gICAgdmFyIGlzID0gdGhpcy5leHRyYUluZGljZXNbMF07XG4gICAgaXMubG9va3VwKHF1ZXJ5LCBmdW5jdGlvbihyZXN1bHQsIGVycikge1xuICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogcmVzdWx0LCBlcnJvcjogZXJyfSk7XG4gICAgfSk7XG59XG5cbn0pLmNhbGwodGhpcyx0eXBlb2Ygc2VsZiAhPT0gXCJ1bmRlZmluZWRcIiA/IHNlbGYgOiB0eXBlb2Ygd2luZG93ICE9PSBcInVuZGVmaW5lZFwiID8gd2luZG93IDoge30pIiwiLyogLSotIG1vZGU6IGphdmFzY3JpcHQ7IGMtYmFzaWMtb2Zmc2V0OiA0OyBpbmRlbnQtdGFicy1tb2RlOiBuaWwgLSotICovXG5cbi8vIFxuLy8gRGFsbGlhbmNlIEdlbm9tZSBFeHBsb3JlclxuLy8gKGMpIFRob21hcyBEb3duIDIwMDYtMjAxMVxuLy9cbi8vIGxoM3V0aWxzLmpzOiBjb21tb24gc3VwcG9ydCBmb3IgbGgzJ3MgZmlsZSBmb3JtYXRzXG4vL1xuXG5pZiAodHlwZW9mKHJlcXVpcmUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIHZhciBqc3psaWIgPSByZXF1aXJlKCdqc3psaWInKTtcbiAgICB2YXIganN6bGliX2luZmxhdGVfYnVmZmVyID0ganN6bGliLmluZmxhdGVCdWZmZXI7XG4gICAgdmFyIGFycmF5Q29weSA9IGpzemxpYi5hcnJheUNvcHk7XG59XG5cbmZ1bmN0aW9uIFZvYihiLCBvKSB7XG4gICAgdGhpcy5ibG9jayA9IGI7XG4gICAgdGhpcy5vZmZzZXQgPSBvO1xufVxuXG5Wb2IucHJvdG90eXBlLnRvU3RyaW5nID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuICcnICsgdGhpcy5ibG9jayArICc6JyArIHRoaXMub2Zmc2V0O1xufVxuXG5mdW5jdGlvbiByZWFkVm9iKGJhLCBvZmZzZXQpIHtcbiAgICB2YXIgYmxvY2sgPSAoKGJhW29mZnNldCs2XSAmIDB4ZmYpICogMHgxMDAwMDAwMDApICsgKChiYVtvZmZzZXQrNV0gJiAweGZmKSAqIDB4MTAwMDAwMCkgKyAoKGJhW29mZnNldCs0XSAmIDB4ZmYpICogMHgxMDAwMCkgKyAoKGJhW29mZnNldCszXSAmIDB4ZmYpICogMHgxMDApICsgKChiYVtvZmZzZXQrMl0gJiAweGZmKSk7XG4gICAgdmFyIGJpbnQgPSAoYmFbb2Zmc2V0KzFdIDw8IDgpIHwgKGJhW29mZnNldF0pO1xuICAgIGlmIChibG9jayA9PSAwICYmIGJpbnQgPT0gMCkge1xuICAgICAgICByZXR1cm4gbnVsbDsgIC8vIFNob3VsZCBvbmx5IGhhcHBlbiBpbiB0aGUgbGluZWFyIGluZGV4P1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHJldHVybiBuZXcgVm9iKGJsb2NrLCBiaW50KTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIHVuYmd6ZihkYXRhLCBsaW0pIHtcbiAgICBsaW0gPSBNYXRoLm1pbihsaW0gfHwgMSwgZGF0YS5ieXRlTGVuZ3RoIC0gNTApO1xuICAgIHZhciBvQmxvY2tMaXN0ID0gW107XG4gICAgdmFyIHB0ciA9IFswXTtcbiAgICB2YXIgdG90YWxTaXplID0gMDtcblxuICAgIHdoaWxlIChwdHJbMF0gPCBsaW0pIHtcbiAgICAgICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoZGF0YSwgcHRyWzBdLCAxMik7IC8vIEZJWE1FIGlzIHRoaXMgZW5vdWdoIGZvciBhbGwgY3JlZGlibGUgQkdaRiBibG9jayBoZWFkZXJzP1xuICAgICAgICB2YXIgeGxlbiA9IChiYVsxMV0gPDwgOCkgfCAoYmFbMTBdKTtcbiAgICAgICAgLy8gZGxvZygneGxlblsnICsgKHB0clswXSkgKyddPScgKyB4bGVuKTtcbiAgICAgICAgdmFyIHVuYyA9IGpzemxpYl9pbmZsYXRlX2J1ZmZlcihkYXRhLCAxMiArIHhsZW4gKyBwdHJbMF0sIE1hdGgubWluKDY1NTM2LCBkYXRhLmJ5dGVMZW5ndGggLSAxMiAtIHhsZW4gLSBwdHJbMF0pLCBwdHIpO1xuICAgICAgICBwdHJbMF0gKz0gODtcbiAgICAgICAgdG90YWxTaXplICs9IHVuYy5ieXRlTGVuZ3RoO1xuICAgICAgICBvQmxvY2tMaXN0LnB1c2godW5jKTtcbiAgICB9XG5cbiAgICBpZiAob0Jsb2NrTGlzdC5sZW5ndGggPT0gMSkge1xuICAgICAgICByZXR1cm4gb0Jsb2NrTGlzdFswXTtcbiAgICB9IGVsc2Uge1xuICAgICAgICB2YXIgb3V0ID0gbmV3IFVpbnQ4QXJyYXkodG90YWxTaXplKTtcbiAgICAgICAgdmFyIGN1cnNvciA9IDA7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgb0Jsb2NrTGlzdC5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgdmFyIGIgPSBuZXcgVWludDhBcnJheShvQmxvY2tMaXN0W2ldKTtcbiAgICAgICAgICAgIGFycmF5Q29weShiLCAwLCBvdXQsIGN1cnNvciwgYi5sZW5ndGgpO1xuICAgICAgICAgICAgY3Vyc29yICs9IGIubGVuZ3RoO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBvdXQuYnVmZmVyO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gQ2h1bmsobWludiwgbWF4dikge1xuICAgIHRoaXMubWludiA9IG1pbnY7IHRoaXMubWF4diA9IG1heHY7XG59XG5cblxuLy9cbi8vIEJpbm5pbmcgKHRyYW5zbGl0ZXJhdGVkIGZyb20gU0FNMS4zIHNwZWMpXG4vL1xuXG4vKiBjYWxjdWxhdGUgYmluIGdpdmVuIGFuIGFsaWdubWVudCBjb3ZlcmluZyBbYmVnLGVuZCkgKHplcm8tYmFzZWQsIGhhbGYtY2xvc2UtaGFsZi1vcGVuKSAqL1xuZnVuY3Rpb24gcmVnMmJpbihiZWcsIGVuZClcbntcbiAgICAtLWVuZDtcbiAgICBpZiAoYmVnPj4xNCA9PSBlbmQ+PjE0KSByZXR1cm4gKCgxPDwxNSktMSkvNyArIChiZWc+PjE0KTtcbiAgICBpZiAoYmVnPj4xNyA9PSBlbmQ+PjE3KSByZXR1cm4gKCgxPDwxMiktMSkvNyArIChiZWc+PjE3KTtcbiAgICBpZiAoYmVnPj4yMCA9PSBlbmQ+PjIwKSByZXR1cm4gKCgxPDw5KS0xKS83ICsgKGJlZz4+MjApO1xuICAgIGlmIChiZWc+PjIzID09IGVuZD4+MjMpIHJldHVybiAoKDE8PDYpLTEpLzcgKyAoYmVnPj4yMyk7XG4gICAgaWYgKGJlZz4+MjYgPT0gZW5kPj4yNikgcmV0dXJuICgoMTw8MyktMSkvNyArIChiZWc+PjI2KTtcbiAgICByZXR1cm4gMDtcbn1cblxuLyogY2FsY3VsYXRlIHRoZSBsaXN0IG9mIGJpbnMgdGhhdCBtYXkgb3ZlcmxhcCB3aXRoIHJlZ2lvbiBbYmVnLGVuZCkgKHplcm8tYmFzZWQpICovXG52YXIgTUFYX0JJTiA9ICgoKDE8PDE4KS0xKS83KTtcbmZ1bmN0aW9uIHJlZzJiaW5zKGJlZywgZW5kKSBcbntcbiAgICB2YXIgaSA9IDAsIGssIGxpc3QgPSBbXTtcbiAgICAtLWVuZDtcbiAgICBsaXN0LnB1c2goMCk7XG4gICAgZm9yIChrID0gMSArIChiZWc+PjI2KTsgayA8PSAxICsgKGVuZD4+MjYpOyArK2spIGxpc3QucHVzaChrKTtcbiAgICBmb3IgKGsgPSA5ICsgKGJlZz4+MjMpOyBrIDw9IDkgKyAoZW5kPj4yMyk7ICsraykgbGlzdC5wdXNoKGspO1xuICAgIGZvciAoayA9IDczICsgKGJlZz4+MjApOyBrIDw9IDczICsgKGVuZD4+MjApOyArK2spIGxpc3QucHVzaChrKTtcbiAgICBmb3IgKGsgPSA1ODUgKyAoYmVnPj4xNyk7IGsgPD0gNTg1ICsgKGVuZD4+MTcpOyArK2spIGxpc3QucHVzaChrKTtcbiAgICBmb3IgKGsgPSA0NjgxICsgKGJlZz4+MTQpOyBrIDw9IDQ2ODEgKyAoZW5kPj4xNCk7ICsraykgbGlzdC5wdXNoKGspO1xuICAgIHJldHVybiBsaXN0O1xufVxuXG5pZiAodHlwZW9mKG1vZHVsZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgICAgIHVuYmd6ZjogdW5iZ3pmLFxuICAgICAgICByZWFkVm9iOiByZWFkVm9iLFxuICAgICAgICByZWcyYmluOiByZWcyYmluLFxuICAgICAgICByZWcyYmluczogcmVnMmJpbnMsXG4gICAgICAgIENodW5rOiBDaHVua1xuICAgIH07XG59IiwiLypcclxuICogQSBKYXZhU2NyaXB0IGltcGxlbWVudGF0aW9uIG9mIHRoZSBTZWN1cmUgSGFzaCBBbGdvcml0aG0sIFNIQS0xLCBhcyBkZWZpbmVkXHJcbiAqIGluIEZJUFMgMTgwLTFcclxuICogVmVyc2lvbiAyLjIgQ29weXJpZ2h0IFBhdWwgSm9obnN0b24gMjAwMCAtIDIwMDkuXHJcbiAqIE90aGVyIGNvbnRyaWJ1dG9yczogR3JlZyBIb2x0LCBBbmRyZXcgS2VwZXJ0LCBZZG5hciwgTG9zdGluZXRcclxuICogRGlzdHJpYnV0ZWQgdW5kZXIgdGhlIEJTRCBMaWNlbnNlXHJcbiAqIFNlZSBodHRwOi8vcGFqaG9tZS5vcmcudWsvY3J5cHQvbWQ1IGZvciBkZXRhaWxzLlxyXG4gKi9cclxuXHJcbiBcInVzZSBzdHJpY3RcIjtcclxuXHJcbi8qXHJcbiAqIENvbmZpZ3VyYWJsZSB2YXJpYWJsZXMuIFlvdSBtYXkgbmVlZCB0byB0d2VhayB0aGVzZSB0byBiZSBjb21wYXRpYmxlIHdpdGhcclxuICogdGhlIHNlcnZlci1zaWRlLCBidXQgdGhlIGRlZmF1bHRzIHdvcmsgaW4gbW9zdCBjYXNlcy5cclxuICovXHJcbnZhciBoZXhjYXNlID0gMDsgIC8qIGhleCBvdXRwdXQgZm9ybWF0LiAwIC0gbG93ZXJjYXNlOyAxIC0gdXBwZXJjYXNlICAgICAgICAqL1xyXG52YXIgYjY0cGFkICA9IFwiXCI7IC8qIGJhc2UtNjQgcGFkIGNoYXJhY3Rlci4gXCI9XCIgZm9yIHN0cmljdCBSRkMgY29tcGxpYW5jZSAgICovXHJcblxyXG4vKlxyXG4gKiBUaGVzZSBhcmUgdGhlIGZ1bmN0aW9ucyB5b3UnbGwgdXN1YWxseSB3YW50IHRvIGNhbGxcclxuICogVGhleSB0YWtlIHN0cmluZyBhcmd1bWVudHMgYW5kIHJldHVybiBlaXRoZXIgaGV4IG9yIGJhc2UtNjQgZW5jb2RlZCBzdHJpbmdzXHJcbiAqL1xyXG5mdW5jdGlvbiBoZXhfc2hhMShzKSAgICB7IHJldHVybiByc3RyMmhleChyc3RyX3NoYTEoc3RyMnJzdHJfdXRmOChzKSkpOyB9XHJcbmZ1bmN0aW9uIGI2NF9zaGExKHMpICAgIHsgcmV0dXJuIHJzdHIyYjY0KHJzdHJfc2hhMShzdHIycnN0cl91dGY4KHMpKSk7IH1cclxuZnVuY3Rpb24gYW55X3NoYTEocywgZSkgeyByZXR1cm4gcnN0cjJhbnkocnN0cl9zaGExKHN0cjJyc3RyX3V0ZjgocykpLCBlKTsgfVxyXG5mdW5jdGlvbiBoZXhfaG1hY19zaGExKGssIGQpXHJcbiAgeyByZXR1cm4gcnN0cjJoZXgocnN0cl9obWFjX3NoYTEoc3RyMnJzdHJfdXRmOChrKSwgc3RyMnJzdHJfdXRmOChkKSkpOyB9XHJcbmZ1bmN0aW9uIGI2NF9obWFjX3NoYTEoaywgZClcclxuICB7IHJldHVybiByc3RyMmI2NChyc3RyX2htYWNfc2hhMShzdHIycnN0cl91dGY4KGspLCBzdHIycnN0cl91dGY4KGQpKSk7IH1cclxuZnVuY3Rpb24gYW55X2htYWNfc2hhMShrLCBkLCBlKVxyXG4gIHsgcmV0dXJuIHJzdHIyYW55KHJzdHJfaG1hY19zaGExKHN0cjJyc3RyX3V0ZjgoayksIHN0cjJyc3RyX3V0ZjgoZCkpLCBlKTsgfVxyXG5cclxuLypcclxuICogUGVyZm9ybSBhIHNpbXBsZSBzZWxmLXRlc3QgdG8gc2VlIGlmIHRoZSBWTSBpcyB3b3JraW5nXHJcbiAqL1xyXG5mdW5jdGlvbiBzaGExX3ZtX3Rlc3QoKVxyXG57XHJcbiAgcmV0dXJuIGhleF9zaGExKFwiYWJjXCIpLnRvTG93ZXJDYXNlKCkgPT0gXCJhOTk5M2UzNjQ3MDY4MTZhYmEzZTI1NzE3ODUwYzI2YzljZDBkODlkXCI7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENhbGN1bGF0ZSB0aGUgU0hBMSBvZiBhIHJhdyBzdHJpbmdcclxuICovXHJcbmZ1bmN0aW9uIHJzdHJfc2hhMShzKVxyXG57XHJcbiAgcmV0dXJuIGJpbmIycnN0cihiaW5iX3NoYTEocnN0cjJiaW5iKHMpLCBzLmxlbmd0aCAqIDgpKTtcclxufVxyXG5cclxuLypcclxuICogQ2FsY3VsYXRlIHRoZSBITUFDLVNIQTEgb2YgYSBrZXkgYW5kIHNvbWUgZGF0YSAocmF3IHN0cmluZ3MpXHJcbiAqL1xyXG5mdW5jdGlvbiByc3RyX2htYWNfc2hhMShrZXksIGRhdGEpXHJcbntcclxuICB2YXIgYmtleSA9IHJzdHIyYmluYihrZXkpO1xyXG4gIGlmKGJrZXkubGVuZ3RoID4gMTYpIGJrZXkgPSBiaW5iX3NoYTEoYmtleSwga2V5Lmxlbmd0aCAqIDgpO1xyXG5cclxuICB2YXIgaXBhZCA9IEFycmF5KDE2KSwgb3BhZCA9IEFycmF5KDE2KTtcclxuICBmb3IodmFyIGkgPSAwOyBpIDwgMTY7IGkrKylcclxuICB7XHJcbiAgICBpcGFkW2ldID0gYmtleVtpXSBeIDB4MzYzNjM2MzY7XHJcbiAgICBvcGFkW2ldID0gYmtleVtpXSBeIDB4NUM1QzVDNUM7XHJcbiAgfVxyXG5cclxuICB2YXIgaGFzaCA9IGJpbmJfc2hhMShpcGFkLmNvbmNhdChyc3RyMmJpbmIoZGF0YSkpLCA1MTIgKyBkYXRhLmxlbmd0aCAqIDgpO1xyXG4gIHJldHVybiBiaW5iMnJzdHIoYmluYl9zaGExKG9wYWQuY29uY2F0KGhhc2gpLCA1MTIgKyAxNjApKTtcclxufVxyXG5cclxuLypcclxuICogQ29udmVydCBhIHJhdyBzdHJpbmcgdG8gYSBoZXggc3RyaW5nXHJcbiAqL1xyXG5mdW5jdGlvbiByc3RyMmhleChpbnB1dClcclxue1xyXG4gIC8vIHRyeSB7IGhleGNhc2UgfSBjYXRjaChlKSB7IGhleGNhc2U9MDsgfVxyXG4gIHZhciBoZXhfdGFiID0gaGV4Y2FzZSA/IFwiMDEyMzQ1Njc4OUFCQ0RFRlwiIDogXCIwMTIzNDU2Nzg5YWJjZGVmXCI7XHJcbiAgdmFyIG91dHB1dCA9IFwiXCI7XHJcbiAgdmFyIHg7XHJcbiAgZm9yKHZhciBpID0gMDsgaSA8IGlucHV0Lmxlbmd0aDsgaSsrKVxyXG4gIHtcclxuICAgIHggPSBpbnB1dC5jaGFyQ29kZUF0KGkpO1xyXG4gICAgb3V0cHV0ICs9IGhleF90YWIuY2hhckF0KCh4ID4+PiA0KSAmIDB4MEYpXHJcbiAgICAgICAgICAgKyAgaGV4X3RhYi5jaGFyQXQoIHggICAgICAgICYgMHgwRik7XHJcbiAgfVxyXG4gIHJldHVybiBvdXRwdXQ7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENvbnZlcnQgYSByYXcgc3RyaW5nIHRvIGEgYmFzZS02NCBzdHJpbmdcclxuICovXHJcbmZ1bmN0aW9uIHJzdHIyYjY0KGlucHV0KVxyXG57XHJcbiAgLy8gdHJ5IHsgYjY0cGFkIH0gY2F0Y2goZSkgeyBiNjRwYWQ9Jyc7IH1cclxuICB2YXIgdGFiID0gXCJBQkNERUZHSElKS0xNTk9QUVJTVFVWV1hZWmFiY2RlZmdoaWprbG1ub3BxcnN0dXZ3eHl6MDEyMzQ1Njc4OSsvXCI7XHJcbiAgdmFyIG91dHB1dCA9IFwiXCI7XHJcbiAgdmFyIGxlbiA9IGlucHV0Lmxlbmd0aDtcclxuICBmb3IodmFyIGkgPSAwOyBpIDwgbGVuOyBpICs9IDMpXHJcbiAge1xyXG4gICAgdmFyIHRyaXBsZXQgPSAoaW5wdXQuY2hhckNvZGVBdChpKSA8PCAxNilcclxuICAgICAgICAgICAgICAgIHwgKGkgKyAxIDwgbGVuID8gaW5wdXQuY2hhckNvZGVBdChpKzEpIDw8IDggOiAwKVxyXG4gICAgICAgICAgICAgICAgfCAoaSArIDIgPCBsZW4gPyBpbnB1dC5jaGFyQ29kZUF0KGkrMikgICAgICA6IDApO1xyXG4gICAgZm9yKHZhciBqID0gMDsgaiA8IDQ7IGorKylcclxuICAgIHtcclxuICAgICAgaWYoaSAqIDggKyBqICogNiA+IGlucHV0Lmxlbmd0aCAqIDgpIG91dHB1dCArPSBiNjRwYWQ7XHJcbiAgICAgIGVsc2Ugb3V0cHV0ICs9IHRhYi5jaGFyQXQoKHRyaXBsZXQgPj4+IDYqKDMtaikpICYgMHgzRik7XHJcbiAgICB9XHJcbiAgfVxyXG4gIHJldHVybiBvdXRwdXQ7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENvbnZlcnQgYSByYXcgc3RyaW5nIHRvIGFuIGFyYml0cmFyeSBzdHJpbmcgZW5jb2RpbmdcclxuICovXHJcbmZ1bmN0aW9uIHJzdHIyYW55KGlucHV0LCBlbmNvZGluZylcclxue1xyXG4gIHZhciBkaXZpc29yID0gZW5jb2RpbmcubGVuZ3RoO1xyXG4gIHZhciByZW1haW5kZXJzID0gQXJyYXkoKTtcclxuICB2YXIgaSwgcSwgeCwgcXVvdGllbnQ7XHJcblxyXG4gIC8qIENvbnZlcnQgdG8gYW4gYXJyYXkgb2YgMTYtYml0IGJpZy1lbmRpYW4gdmFsdWVzLCBmb3JtaW5nIHRoZSBkaXZpZGVuZCAqL1xyXG4gIHZhciBkaXZpZGVuZCA9IEFycmF5KE1hdGguY2VpbChpbnB1dC5sZW5ndGggLyAyKSk7XHJcbiAgZm9yKGkgPSAwOyBpIDwgZGl2aWRlbmQubGVuZ3RoOyBpKyspXHJcbiAge1xyXG4gICAgZGl2aWRlbmRbaV0gPSAoaW5wdXQuY2hhckNvZGVBdChpICogMikgPDwgOCkgfCBpbnB1dC5jaGFyQ29kZUF0KGkgKiAyICsgMSk7XHJcbiAgfVxyXG5cclxuICAvKlxyXG4gICAqIFJlcGVhdGVkbHkgcGVyZm9ybSBhIGxvbmcgZGl2aXNpb24uIFRoZSBiaW5hcnkgYXJyYXkgZm9ybXMgdGhlIGRpdmlkZW5kLFxyXG4gICAqIHRoZSBsZW5ndGggb2YgdGhlIGVuY29kaW5nIGlzIHRoZSBkaXZpc29yLiBPbmNlIGNvbXB1dGVkLCB0aGUgcXVvdGllbnRcclxuICAgKiBmb3JtcyB0aGUgZGl2aWRlbmQgZm9yIHRoZSBuZXh0IHN0ZXAuIFdlIHN0b3Agd2hlbiB0aGUgZGl2aWRlbmQgaXMgemVyby5cclxuICAgKiBBbGwgcmVtYWluZGVycyBhcmUgc3RvcmVkIGZvciBsYXRlciB1c2UuXHJcbiAgICovXHJcbiAgd2hpbGUoZGl2aWRlbmQubGVuZ3RoID4gMClcclxuICB7XHJcbiAgICBxdW90aWVudCA9IEFycmF5KCk7XHJcbiAgICB4ID0gMDtcclxuICAgIGZvcihpID0gMDsgaSA8IGRpdmlkZW5kLmxlbmd0aDsgaSsrKVxyXG4gICAge1xyXG4gICAgICB4ID0gKHggPDwgMTYpICsgZGl2aWRlbmRbaV07XHJcbiAgICAgIHEgPSBNYXRoLmZsb29yKHggLyBkaXZpc29yKTtcclxuICAgICAgeCAtPSBxICogZGl2aXNvcjtcclxuICAgICAgaWYocXVvdGllbnQubGVuZ3RoID4gMCB8fCBxID4gMClcclxuICAgICAgICBxdW90aWVudFtxdW90aWVudC5sZW5ndGhdID0gcTtcclxuICAgIH1cclxuICAgIHJlbWFpbmRlcnNbcmVtYWluZGVycy5sZW5ndGhdID0geDtcclxuICAgIGRpdmlkZW5kID0gcXVvdGllbnQ7XHJcbiAgfVxyXG5cclxuICAvKiBDb252ZXJ0IHRoZSByZW1haW5kZXJzIHRvIHRoZSBvdXRwdXQgc3RyaW5nICovXHJcbiAgdmFyIG91dHB1dCA9IFwiXCI7XHJcbiAgZm9yKGkgPSByZW1haW5kZXJzLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKVxyXG4gICAgb3V0cHV0ICs9IGVuY29kaW5nLmNoYXJBdChyZW1haW5kZXJzW2ldKTtcclxuXHJcbiAgLyogQXBwZW5kIGxlYWRpbmcgemVybyBlcXVpdmFsZW50cyAqL1xyXG4gIHZhciBmdWxsX2xlbmd0aCA9IE1hdGguY2VpbChpbnB1dC5sZW5ndGggKiA4IC9cclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgKE1hdGgubG9nKGVuY29kaW5nLmxlbmd0aCkgLyBNYXRoLmxvZygyKSkpXHJcbiAgZm9yKGkgPSBvdXRwdXQubGVuZ3RoOyBpIDwgZnVsbF9sZW5ndGg7IGkrKylcclxuICAgIG91dHB1dCA9IGVuY29kaW5nWzBdICsgb3V0cHV0O1xyXG5cclxuICByZXR1cm4gb3V0cHV0O1xyXG59XHJcblxyXG4vKlxyXG4gKiBFbmNvZGUgYSBzdHJpbmcgYXMgdXRmLTguXHJcbiAqIEZvciBlZmZpY2llbmN5LCB0aGlzIGFzc3VtZXMgdGhlIGlucHV0IGlzIHZhbGlkIHV0Zi0xNi5cclxuICovXHJcbmZ1bmN0aW9uIHN0cjJyc3RyX3V0ZjgoaW5wdXQpXHJcbntcclxuICB2YXIgb3V0cHV0ID0gXCJcIjtcclxuICB2YXIgaSA9IC0xO1xyXG4gIHZhciB4LCB5O1xyXG5cclxuICB3aGlsZSgrK2kgPCBpbnB1dC5sZW5ndGgpXHJcbiAge1xyXG4gICAgLyogRGVjb2RlIHV0Zi0xNiBzdXJyb2dhdGUgcGFpcnMgKi9cclxuICAgIHggPSBpbnB1dC5jaGFyQ29kZUF0KGkpO1xyXG4gICAgeSA9IGkgKyAxIDwgaW5wdXQubGVuZ3RoID8gaW5wdXQuY2hhckNvZGVBdChpICsgMSkgOiAwO1xyXG4gICAgaWYoMHhEODAwIDw9IHggJiYgeCA8PSAweERCRkYgJiYgMHhEQzAwIDw9IHkgJiYgeSA8PSAweERGRkYpXHJcbiAgICB7XHJcbiAgICAgIHggPSAweDEwMDAwICsgKCh4ICYgMHgwM0ZGKSA8PCAxMCkgKyAoeSAmIDB4MDNGRik7XHJcbiAgICAgIGkrKztcclxuICAgIH1cclxuXHJcbiAgICAvKiBFbmNvZGUgb3V0cHV0IGFzIHV0Zi04ICovXHJcbiAgICBpZih4IDw9IDB4N0YpXHJcbiAgICAgIG91dHB1dCArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKHgpO1xyXG4gICAgZWxzZSBpZih4IDw9IDB4N0ZGKVxyXG4gICAgICBvdXRwdXQgKz0gU3RyaW5nLmZyb21DaGFyQ29kZSgweEMwIHwgKCh4ID4+PiA2ICkgJiAweDFGKSxcclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgMHg4MCB8ICggeCAgICAgICAgICYgMHgzRikpO1xyXG4gICAgZWxzZSBpZih4IDw9IDB4RkZGRilcclxuICAgICAgb3V0cHV0ICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoMHhFMCB8ICgoeCA+Pj4gMTIpICYgMHgwRiksXHJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIDB4ODAgfCAoKHggPj4+IDYgKSAmIDB4M0YpLFxyXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAweDgwIHwgKCB4ICAgICAgICAgJiAweDNGKSk7XHJcbiAgICBlbHNlIGlmKHggPD0gMHgxRkZGRkYpXHJcbiAgICAgIG91dHB1dCArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKDB4RjAgfCAoKHggPj4+IDE4KSAmIDB4MDcpLFxyXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAweDgwIHwgKCh4ID4+PiAxMikgJiAweDNGKSxcclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgMHg4MCB8ICgoeCA+Pj4gNiApICYgMHgzRiksXHJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIDB4ODAgfCAoIHggICAgICAgICAmIDB4M0YpKTtcclxuICB9XHJcbiAgcmV0dXJuIG91dHB1dDtcclxufVxyXG5cclxuLypcclxuICogRW5jb2RlIGEgc3RyaW5nIGFzIHV0Zi0xNlxyXG4gKi9cclxuZnVuY3Rpb24gc3RyMnJzdHJfdXRmMTZsZShpbnB1dClcclxue1xyXG4gIHZhciBvdXRwdXQgPSBcIlwiO1xyXG4gIGZvcih2YXIgaSA9IDA7IGkgPCBpbnB1dC5sZW5ndGg7IGkrKylcclxuICAgIG91dHB1dCArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKCBpbnB1dC5jaGFyQ29kZUF0KGkpICAgICAgICAmIDB4RkYsXHJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAoaW5wdXQuY2hhckNvZGVBdChpKSA+Pj4gOCkgJiAweEZGKTtcclxuICByZXR1cm4gb3V0cHV0O1xyXG59XHJcblxyXG5mdW5jdGlvbiBzdHIycnN0cl91dGYxNmJlKGlucHV0KVxyXG57XHJcbiAgdmFyIG91dHB1dCA9IFwiXCI7XHJcbiAgZm9yKHZhciBpID0gMDsgaSA8IGlucHV0Lmxlbmd0aDsgaSsrKVxyXG4gICAgb3V0cHV0ICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoKGlucHV0LmNoYXJDb2RlQXQoaSkgPj4+IDgpICYgMHhGRixcclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbnB1dC5jaGFyQ29kZUF0KGkpICAgICAgICAmIDB4RkYpO1xyXG4gIHJldHVybiBvdXRwdXQ7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENvbnZlcnQgYSByYXcgc3RyaW5nIHRvIGFuIGFycmF5IG9mIGJpZy1lbmRpYW4gd29yZHNcclxuICogQ2hhcmFjdGVycyA+MjU1IGhhdmUgdGhlaXIgaGlnaC1ieXRlIHNpbGVudGx5IGlnbm9yZWQuXHJcbiAqL1xyXG5mdW5jdGlvbiByc3RyMmJpbmIoaW5wdXQpXHJcbntcclxuICB2YXIgb3V0cHV0ID0gQXJyYXkoaW5wdXQubGVuZ3RoID4+IDIpO1xyXG4gIGZvcih2YXIgaSA9IDA7IGkgPCBvdXRwdXQubGVuZ3RoOyBpKyspXHJcbiAgICBvdXRwdXRbaV0gPSAwO1xyXG4gIGZvcih2YXIgaSA9IDA7IGkgPCBpbnB1dC5sZW5ndGggKiA4OyBpICs9IDgpXHJcbiAgICBvdXRwdXRbaT4+NV0gfD0gKGlucHV0LmNoYXJDb2RlQXQoaSAvIDgpICYgMHhGRikgPDwgKDI0IC0gaSAlIDMyKTtcclxuICByZXR1cm4gb3V0cHV0O1xyXG59XHJcblxyXG4vKlxyXG4gKiBDb252ZXJ0IGFuIGFycmF5IG9mIGJpZy1lbmRpYW4gd29yZHMgdG8gYSBzdHJpbmdcclxuICovXHJcbmZ1bmN0aW9uIGJpbmIycnN0cihpbnB1dClcclxue1xyXG4gIHZhciBvdXRwdXQgPSBcIlwiO1xyXG4gIGZvcih2YXIgaSA9IDA7IGkgPCBpbnB1dC5sZW5ndGggKiAzMjsgaSArPSA4KVxyXG4gICAgb3V0cHV0ICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoKGlucHV0W2k+PjVdID4+PiAoMjQgLSBpICUgMzIpKSAmIDB4RkYpO1xyXG4gIHJldHVybiBvdXRwdXQ7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENhbGN1bGF0ZSB0aGUgU0hBLTEgb2YgYW4gYXJyYXkgb2YgYmlnLWVuZGlhbiB3b3JkcywgYW5kIGEgYml0IGxlbmd0aFxyXG4gKi9cclxuZnVuY3Rpb24gYmluYl9zaGExKHgsIGxlbilcclxue1xyXG4gIC8qIGFwcGVuZCBwYWRkaW5nICovXHJcbiAgeFtsZW4gPj4gNV0gfD0gMHg4MCA8PCAoMjQgLSBsZW4gJSAzMik7XHJcbiAgeFsoKGxlbiArIDY0ID4+IDkpIDw8IDQpICsgMTVdID0gbGVuO1xyXG5cclxuICB2YXIgdyA9IEFycmF5KDgwKTtcclxuICB2YXIgYSA9ICAxNzMyNTg0MTkzO1xyXG4gIHZhciBiID0gLTI3MTczMzg3OTtcclxuICB2YXIgYyA9IC0xNzMyNTg0MTk0O1xyXG4gIHZhciBkID0gIDI3MTczMzg3ODtcclxuICB2YXIgZSA9IC0xMDA5NTg5Nzc2O1xyXG5cclxuICBmb3IodmFyIGkgPSAwOyBpIDwgeC5sZW5ndGg7IGkgKz0gMTYpXHJcbiAge1xyXG4gICAgdmFyIG9sZGEgPSBhO1xyXG4gICAgdmFyIG9sZGIgPSBiO1xyXG4gICAgdmFyIG9sZGMgPSBjO1xyXG4gICAgdmFyIG9sZGQgPSBkO1xyXG4gICAgdmFyIG9sZGUgPSBlO1xyXG5cclxuICAgIGZvcih2YXIgaiA9IDA7IGogPCA4MDsgaisrKVxyXG4gICAge1xyXG4gICAgICBpZihqIDwgMTYpIHdbal0gPSB4W2kgKyBqXTtcclxuICAgICAgZWxzZSB3W2pdID0gYml0X3JvbCh3W2otM10gXiB3W2otOF0gXiB3W2otMTRdIF4gd1tqLTE2XSwgMSk7XHJcbiAgICAgIHZhciB0ID0gc2FmZV9hZGQoc2FmZV9hZGQoYml0X3JvbChhLCA1KSwgc2hhMV9mdChqLCBiLCBjLCBkKSksXHJcbiAgICAgICAgICAgICAgICAgICAgICAgc2FmZV9hZGQoc2FmZV9hZGQoZSwgd1tqXSksIHNoYTFfa3QoaikpKTtcclxuICAgICAgZSA9IGQ7XHJcbiAgICAgIGQgPSBjO1xyXG4gICAgICBjID0gYml0X3JvbChiLCAzMCk7XHJcbiAgICAgIGIgPSBhO1xyXG4gICAgICBhID0gdDtcclxuICAgIH1cclxuXHJcbiAgICBhID0gc2FmZV9hZGQoYSwgb2xkYSk7XHJcbiAgICBiID0gc2FmZV9hZGQoYiwgb2xkYik7XHJcbiAgICBjID0gc2FmZV9hZGQoYywgb2xkYyk7XHJcbiAgICBkID0gc2FmZV9hZGQoZCwgb2xkZCk7XHJcbiAgICBlID0gc2FmZV9hZGQoZSwgb2xkZSk7XHJcbiAgfVxyXG4gIHJldHVybiBBcnJheShhLCBiLCBjLCBkLCBlKTtcclxuXHJcbn1cclxuXHJcbi8qXHJcbiAqIFBlcmZvcm0gdGhlIGFwcHJvcHJpYXRlIHRyaXBsZXQgY29tYmluYXRpb24gZnVuY3Rpb24gZm9yIHRoZSBjdXJyZW50XHJcbiAqIGl0ZXJhdGlvblxyXG4gKi9cclxuZnVuY3Rpb24gc2hhMV9mdCh0LCBiLCBjLCBkKVxyXG57XHJcbiAgaWYodCA8IDIwKSByZXR1cm4gKGIgJiBjKSB8ICgofmIpICYgZCk7XHJcbiAgaWYodCA8IDQwKSByZXR1cm4gYiBeIGMgXiBkO1xyXG4gIGlmKHQgPCA2MCkgcmV0dXJuIChiICYgYykgfCAoYiAmIGQpIHwgKGMgJiBkKTtcclxuICByZXR1cm4gYiBeIGMgXiBkO1xyXG59XHJcblxyXG4vKlxyXG4gKiBEZXRlcm1pbmUgdGhlIGFwcHJvcHJpYXRlIGFkZGl0aXZlIGNvbnN0YW50IGZvciB0aGUgY3VycmVudCBpdGVyYXRpb25cclxuICovXHJcbmZ1bmN0aW9uIHNoYTFfa3QodClcclxue1xyXG4gIHJldHVybiAodCA8IDIwKSA/ICAxNTE4NTAwMjQ5IDogKHQgPCA0MCkgPyAgMTg1OTc3NTM5MyA6XHJcbiAgICAgICAgICh0IDwgNjApID8gLTE4OTQwMDc1ODggOiAtODk5NDk3NTE0O1xyXG59XHJcblxyXG4vKlxyXG4gKiBBZGQgaW50ZWdlcnMsIHdyYXBwaW5nIGF0IDJeMzIuIFRoaXMgdXNlcyAxNi1iaXQgb3BlcmF0aW9ucyBpbnRlcm5hbGx5XHJcbiAqIHRvIHdvcmsgYXJvdW5kIGJ1Z3MgaW4gc29tZSBKUyBpbnRlcnByZXRlcnMuXHJcbiAqL1xyXG5mdW5jdGlvbiBzYWZlX2FkZCh4LCB5KVxyXG57XHJcbiAgdmFyIGxzdyA9ICh4ICYgMHhGRkZGKSArICh5ICYgMHhGRkZGKTtcclxuICB2YXIgbXN3ID0gKHggPj4gMTYpICsgKHkgPj4gMTYpICsgKGxzdyA+PiAxNik7XHJcbiAgcmV0dXJuIChtc3cgPDwgMTYpIHwgKGxzdyAmIDB4RkZGRik7XHJcbn1cclxuXHJcbi8qXHJcbiAqIEJpdHdpc2Ugcm90YXRlIGEgMzItYml0IG51bWJlciB0byB0aGUgbGVmdC5cclxuICovXHJcbmZ1bmN0aW9uIGJpdF9yb2wobnVtLCBjbnQpXHJcbntcclxuICByZXR1cm4gKG51bSA8PCBjbnQpIHwgKG51bSA+Pj4gKDMyIC0gY250KSk7XHJcbn1cclxuXHJcbmlmICh0eXBlb2YobW9kdWxlKSAhPT0gJ3VuZGVmaW5lZCcpIHtcclxuICBtb2R1bGUuZXhwb3J0cyA9IHtcclxuICAgIGI2NF9zaGExOiBiNjRfc2hhMSxcclxuICAgIGhleF9zaGExOiBoZXhfc2hhMVxyXG4gIH1cclxufVxyXG4iLCIvKiAtKi0gbW9kZTogamF2YXNjcmlwdDsgYy1iYXNpYy1vZmZzZXQ6IDQ7IGluZGVudC10YWJzLW1vZGU6IG5pbCAtKi0gKi9cblxuLy8gXG4vLyBEYWxsaWFuY2UgR2Vub21lIEV4cGxvcmVyXG4vLyAoYykgVGhvbWFzIERvd24gMjAwNi0yMDEwXG4vL1xuLy8gc3BhbnMuanM6IEphdmFTY3JpcHQgSW50c2V0L0xvY2F0aW9uIHBvcnQuXG4vL1xuXG5cInVzZSBzdHJpY3RcIjtcblxuXG5mdW5jdGlvbiBSYW5nZShtaW4sIG1heClcbntcbiAgICBpZiAodHlwZW9mKG1pbikgIT0gJ251bWJlcicgfHwgdHlwZW9mKG1heCkgIT0gJ251bWJlcicpXG4gICAgICAgIHRocm93ICdCYWQgcmFuZ2UgJyArIG1pbiArICcsJyArIG1heDtcbiAgICB0aGlzLl9taW4gPSBtaW47XG4gICAgdGhpcy5fbWF4ID0gbWF4O1xufVxuXG5SYW5nZS5wcm90b3R5cGUubWluID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX21pbjtcbn1cblxuUmFuZ2UucHJvdG90eXBlLm1heCA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9tYXg7XG59XG5cblJhbmdlLnByb3RvdHlwZS5jb250YWlucyA9IGZ1bmN0aW9uKHBvcykge1xuICAgIHJldHVybiBwb3MgPj0gdGhpcy5fbWluICYmIHBvcyA8PSB0aGlzLl9tYXg7XG59XG5cblJhbmdlLnByb3RvdHlwZS5pc0NvbnRpZ3VvdXMgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdHJ1ZTtcbn1cblxuUmFuZ2UucHJvdG90eXBlLnJhbmdlcyA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiBbdGhpc107XG59XG5cblJhbmdlLnByb3RvdHlwZS5fcHVzaFJhbmdlcyA9IGZ1bmN0aW9uKHJhbmdlcykge1xuICAgIHJhbmdlcy5wdXNoKHRoaXMpO1xufVxuXG5SYW5nZS5wcm90b3R5cGUudG9TdHJpbmcgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gJ1snICsgdGhpcy5fbWluICsgJy0nICsgdGhpcy5fbWF4ICsgJ10nO1xufVxuXG5mdW5jdGlvbiBfQ29tcG91bmQocmFuZ2VzKSB7XG4gICAgdGhpcy5fcmFuZ2VzID0gcmFuZ2VzO1xuICAgIC8vIGFzc2VydCBzb3J0ZWQ/XG59XG5cbl9Db21wb3VuZC5wcm90b3R5cGUubWluID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX3Jhbmdlc1swXS5taW4oKTtcbn1cblxuX0NvbXBvdW5kLnByb3RvdHlwZS5tYXggPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5fcmFuZ2VzW3RoaXMuX3Jhbmdlcy5sZW5ndGggLSAxXS5tYXgoKTtcbn1cblxuX0NvbXBvdW5kLnByb3RvdHlwZS5jb250YWlucyA9IGZ1bmN0aW9uKHBvcykge1xuICAgIC8vIEZJWE1FIGltcGxlbWVudCBic2VhcmNoIGlmIHdlIHVzZSB0aGlzIG11Y2guXG4gICAgZm9yICh2YXIgcyA9IDA7IHMgPCB0aGlzLl9yYW5nZXMubGVuZ3RoOyArK3MpIHtcbiAgICAgICAgaWYgKHRoaXMuX3Jhbmdlc1tzXS5jb250YWlucyhwb3MpKSB7XG4gICAgICAgICAgICByZXR1cm4gdHJ1ZTtcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gZmFsc2U7XG59XG5cbl9Db21wb3VuZC5wcm90b3R5cGUuaXNDb250aWd1b3VzID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX3Jhbmdlcy5sZW5ndGggPiAxO1xufVxuXG5fQ29tcG91bmQucHJvdG90eXBlLnJhbmdlcyA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9yYW5nZXM7XG59XG5cbl9Db21wb3VuZC5wcm90b3R5cGUuX3B1c2hSYW5nZXMgPSBmdW5jdGlvbihyYW5nZXMpIHtcbiAgICBmb3IgKHZhciByaSA9IDA7IHJpIDwgdGhpcy5fcmFuZ2VzLmxlbmd0aDsgKytyaSlcbiAgICAgICAgcmFuZ2VzLnB1c2godGhpcy5fcmFuZ2VzW3JpXSk7XG59XG5cbl9Db21wb3VuZC5wcm90b3R5cGUudG9TdHJpbmcgPSBmdW5jdGlvbigpIHtcbiAgICB2YXIgcyA9ICcnO1xuICAgIGZvciAodmFyIHIgPSAwOyByIDwgdGhpcy5fcmFuZ2VzLmxlbmd0aDsgKytyKSB7XG4gICAgICAgIGlmIChyPjApIHtcbiAgICAgICAgICAgIHMgPSBzICsgJywnO1xuICAgICAgICB9XG4gICAgICAgIHMgPSBzICsgdGhpcy5fcmFuZ2VzW3JdLnRvU3RyaW5nKCk7XG4gICAgfVxuICAgIHJldHVybiBzO1xufVxuXG5mdW5jdGlvbiB1bmlvbihzMCwgczEpIHtcbiAgICBpZiAoISAoczAgaW5zdGFuY2VvZiBBcnJheSkpIHtcbiAgICAgICAgczAgPSBbczBdO1xuICAgICAgICBpZiAoczEpXG4gICAgICAgICAgICBzMC5wdXNoKHMxKTtcbiAgICB9XG5cbiAgICBpZiAoczAubGVuZ3RoID09IDApXG4gICAgICAgIHJldHVybiBudWxsO1xuICAgIGVsc2UgaWYgKHMwLmxlbmd0aCA9PSAxKVxuICAgICAgICByZXR1cm4gczBbMF07XG5cbiAgICB2YXIgcmFuZ2VzID0gW107XG4gICAgZm9yICh2YXIgc2kgPSAwOyBzaSA8IHMwLmxlbmd0aDsgKytzaSlcbiAgICAgICAgczBbc2ldLl9wdXNoUmFuZ2VzKHJhbmdlcyk7XG4gICAgcmFuZ2VzID0gcmFuZ2VzLnNvcnQoX3JhbmdlT3JkZXIpO1xuXG4gICAgdmFyIG9yYW5nZXMgPSBbXTtcbiAgICB2YXIgY3VycmVudCA9IHJhbmdlc1swXTtcbiAgICBjdXJyZW50ID0gbmV3IFJhbmdlKGN1cnJlbnQuX21pbiwgY3VycmVudC5fbWF4KTsgIC8vIENvcHkgbm93IHNvIHdlIGRvbid0IGhhdmUgdG8gbGF0ZXIuXG5cbiAgICBmb3IgKHZhciBpID0gMTsgaSA8IHJhbmdlcy5sZW5ndGg7ICsraSkge1xuICAgICAgICB2YXIgbnh0ID0gcmFuZ2VzW2ldO1xuICAgICAgICBpZiAobnh0Ll9taW4gPiAoY3VycmVudC5fbWF4ICsgMSkpIHtcbiAgICAgICAgICAgIG9yYW5nZXMucHVzaChjdXJyZW50KTtcbiAgICAgICAgICAgIGN1cnJlbnQgPSBuZXcgUmFuZ2Uobnh0Ll9taW4sIG54dC5fbWF4KTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGlmIChueHQuX21heCA+IGN1cnJlbnQuX21heCkge1xuICAgICAgICAgICAgICAgIGN1cnJlbnQuX21heCA9IG54dC5fbWF4O1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfVxuICAgIG9yYW5nZXMucHVzaChjdXJyZW50KTtcblxuICAgIGlmIChvcmFuZ2VzLmxlbmd0aCA9PSAxKSB7XG4gICAgICAgIHJldHVybiBvcmFuZ2VzWzBdO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHJldHVybiBuZXcgX0NvbXBvdW5kKG9yYW5nZXMpO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gaW50ZXJzZWN0aW9uKHMwLCBzMSkge1xuICAgIHZhciByMCA9IHMwLnJhbmdlcygpO1xuICAgIHZhciByMSA9IHMxLnJhbmdlcygpO1xuICAgIHZhciBsMCA9IHIwLmxlbmd0aCwgbDEgPSByMS5sZW5ndGg7XG4gICAgdmFyIGkwID0gMCwgaTEgPSAwO1xuICAgIHZhciBvciA9IFtdO1xuXG4gICAgd2hpbGUgKGkwIDwgbDAgJiYgaTEgPCBsMSkge1xuICAgICAgICB2YXIgczAgPSByMFtpMF0sIHMxID0gcjFbaTFdO1xuICAgICAgICB2YXIgbGFwTWluID0gTWF0aC5tYXgoczAubWluKCksIHMxLm1pbigpKTtcbiAgICAgICAgdmFyIGxhcE1heCA9IE1hdGgubWluKHMwLm1heCgpLCBzMS5tYXgoKSk7XG4gICAgICAgIGlmIChsYXBNYXggPj0gbGFwTWluKSB7XG4gICAgICAgICAgICBvci5wdXNoKG5ldyBSYW5nZShsYXBNaW4sIGxhcE1heCkpO1xuICAgICAgICB9XG4gICAgICAgIGlmIChzMC5tYXgoKSA+IHMxLm1heCgpKSB7XG4gICAgICAgICAgICArK2kxO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgKytpMDtcbiAgICAgICAgfVxuICAgIH1cbiAgICBcbiAgICBpZiAob3IubGVuZ3RoID09IDApIHtcbiAgICAgICAgcmV0dXJuIG51bGw7IC8vIEZJWE1FXG4gICAgfSBlbHNlIGlmIChvci5sZW5ndGggPT0gMSkge1xuICAgICAgICByZXR1cm4gb3JbMF07XG4gICAgfSBlbHNlIHtcbiAgICAgICAgcmV0dXJuIG5ldyBfQ29tcG91bmQob3IpO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gY292ZXJhZ2Uocykge1xuICAgIHZhciB0b3QgPSAwO1xuICAgIHZhciBybCA9IHMucmFuZ2VzKCk7XG4gICAgZm9yICh2YXIgcmkgPSAwOyByaSA8IHJsLmxlbmd0aDsgKytyaSkge1xuICAgICAgICB2YXIgciA9IHJsW3JpXTtcbiAgICAgICAgdG90ICs9IChyLm1heCgpIC0gci5taW4oKSArIDEpO1xuICAgIH1cbiAgICByZXR1cm4gdG90O1xufVxuXG5cblxuZnVuY3Rpb24gcmFuZ2VPcmRlcihhLCBiKVxue1xuICAgIGlmIChhLm1pbigpIDwgYi5taW4oKSkge1xuICAgICAgICByZXR1cm4gLTE7XG4gICAgfSBlbHNlIGlmIChhLm1pbigpID4gYi5taW4oKSkge1xuICAgICAgICByZXR1cm4gMTtcbiAgICB9IGVsc2UgaWYgKGEubWF4KCkgPCBiLm1heCgpKSB7XG4gICAgICAgIHJldHVybiAtMTtcbiAgICB9IGVsc2UgaWYgKGIubWF4KCkgPiBhLm1heCgpKSB7XG4gICAgICAgIHJldHVybiAxO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHJldHVybiAwO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gX3JhbmdlT3JkZXIoYSwgYilcbntcbiAgICBpZiAoYS5fbWluIDwgYi5fbWluKSB7XG4gICAgICAgIHJldHVybiAtMTtcbiAgICB9IGVsc2UgaWYgKGEuX21pbiA+IGIuX21pbikge1xuICAgICAgICByZXR1cm4gMTtcbiAgICB9IGVsc2UgaWYgKGEuX21heCA8IGIuX21heCkge1xuICAgICAgICByZXR1cm4gLTE7XG4gICAgfSBlbHNlIGlmIChiLl9tYXggPiBhLl9tYXgpIHtcbiAgICAgICAgcmV0dXJuIDE7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgcmV0dXJuIDA7XG4gICAgfVxufVxuXG5pZiAodHlwZW9mKG1vZHVsZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgICAgIFJhbmdlOiBSYW5nZSxcbiAgICAgICAgdW5pb246IHVuaW9uLFxuICAgICAgICBpbnRlcnNlY3Rpb246IGludGVyc2VjdGlvbixcbiAgICAgICAgY292ZXJhZ2U6IGNvdmVyYWdlLFxuICAgICAgICByYW5nZU92ZXI6IHJhbmdlT3JkZXIsXG4gICAgICAgIF9yYW5nZU9yZGVyOiBfcmFuZ2VPcmRlclxuICAgIH1cbn0iLCIvKiAtKi0gbW9kZTogamF2YXNjcmlwdDsgYy1iYXNpYy1vZmZzZXQ6IDQ7IGluZGVudC10YWJzLW1vZGU6IG5pbCAtKi0gKi9cblxuLy8gXG4vLyBEYWxsaWFuY2UgR2Vub21lIEV4cGxvcmVyXG4vLyAoYykgVGhvbWFzIERvd24gMjAwNi0yMDEwXG4vL1xuLy8gdXRpbHMuanM6IG9kZHMsIHNvZHMsIGFuZCBlbmRzLlxuLy9cblxuXCJ1c2Ugc3RyaWN0XCI7XG5cbmlmICh0eXBlb2YocmVxdWlyZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgdmFyIHNoYTEgPSByZXF1aXJlKCcuL3NoYTEnKTtcbiAgICB2YXIgYjY0X3NoYTEgPSBzaGExLmI2NF9zaGExO1xufVxuXG52YXIgTlVNX1JFR0VYUCA9IG5ldyBSZWdFeHAoJ1swLTldKycpO1xuXG5mdW5jdGlvbiBzdHJpbmdUb051bWJlcnNBcnJheShzdHIpIHtcbiAgICB2YXIgbnVtcyA9IG5ldyBBcnJheSgpO1xuICAgIHZhciBtO1xuICAgIHdoaWxlIChtID0gTlVNX1JFR0VYUC5leGVjKHN0cikpIHtcbiAgICAgICAgbnVtcy5wdXNoKG1bMF0pO1xuICAgICAgICBzdHI9c3RyLnN1YnN0cmluZyhtLmluZGV4ICsgKG1bMF0ubGVuZ3RoKSk7XG4gICAgfVxuICAgIHJldHVybiBudW1zO1xufVxuXG52YXIgU1RSSUNUX05VTV9SRUdFWFAgPSBuZXcgUmVnRXhwKCdeWzAtOV0rJCcpO1xuXG5mdW5jdGlvbiBzdHJpbmdUb0ludChzdHIpIHtcbiAgICBzdHIgPSBzdHIucmVwbGFjZShuZXcgUmVnRXhwKCcsJywgJ2cnKSwgJycpO1xuICAgIGlmICghU1RSSUNUX05VTV9SRUdFWFAudGVzdChzdHIpKSB7XG4gICAgICAgIHJldHVybiBudWxsO1xuICAgIH1cbiAgICByZXR1cm4gc3RyfDA7XG59XG5cbmZ1bmN0aW9uIHB1c2huZXcoYSwgdikge1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgYS5sZW5ndGg7ICsraSkge1xuICAgICAgICBpZiAoYVtpXSA9PSB2KSB7XG4gICAgICAgICAgICByZXR1cm47XG4gICAgICAgIH1cbiAgICB9XG4gICAgYS5wdXNoKHYpO1xufVxuXG5mdW5jdGlvbiBwdXNobyhvYmosIGssIHYpIHtcbiAgICBpZiAob2JqW2tdKSB7XG4gICAgICAgIG9ialtrXS5wdXNoKHYpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIG9ialtrXSA9IFt2XTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIHB1c2huZXdvKG9iaiwgaywgdikge1xuICAgIHZhciBhID0gb2JqW2tdO1xuICAgIGlmIChhKSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgYS5sZW5ndGg7ICsraSkgeyAgICAvLyBpbmRleE9mIHJlcXVpcmVzIEpTMTYgOi0oLlxuICAgICAgICAgICAgaWYgKGFbaV0gPT0gdikge1xuICAgICAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBhLnB1c2godik7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgb2JqW2tdID0gW3ZdO1xuICAgIH1cbn1cblxuXG5mdW5jdGlvbiBwaWNrKGEsIGIsIGMsIGQpXG57XG4gICAgaWYgKGEpIHtcbiAgICAgICAgcmV0dXJuIGE7XG4gICAgfSBlbHNlIGlmIChiKSB7XG4gICAgICAgIHJldHVybiBiO1xuICAgIH0gZWxzZSBpZiAoYykge1xuICAgICAgICByZXR1cm4gYztcbiAgICB9IGVsc2UgaWYgKGQpIHtcbiAgICAgICAgcmV0dXJuIGQ7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBwdXNobmV3KGwsIG8pXG57XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBsLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIGlmIChsW2ldID09IG8pIHtcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuICAgIH1cbiAgICBsLnB1c2gobyk7XG59XG5cblxuXG5mdW5jdGlvbiBhcnJheUluZGV4T2YoYSwgeCkge1xuICAgIGlmICghYSkge1xuICAgICAgICByZXR1cm4gLTE7XG4gICAgfVxuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBhLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIGlmIChhW2ldID09PSB4KSB7XG4gICAgICAgICAgICByZXR1cm4gaTtcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gLTE7XG59XG5cbmZ1bmN0aW9uIGFycmF5UmVtb3ZlKGEsIHgpIHtcbiAgICB2YXIgaSA9IGFycmF5SW5kZXhPZihhLCB4KTtcbiAgICBpZiAoaSA+PSAwKSB7XG4gICAgICAgIGEuc3BsaWNlKGksIDEpO1xuICAgICAgICByZXR1cm4gdHJ1ZTtcbiAgICB9XG4gICAgcmV0dXJuIGZhbHNlO1xufVxuXG4vL1xuLy8gRE9NIHV0aWxpdGllc1xuLy9cblxuXG5mdW5jdGlvbiBtYWtlRWxlbWVudCh0YWcsIGNoaWxkcmVuLCBhdHRyaWJzLCBzdHlsZXMpXG57XG4gICAgdmFyIGVsZSA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnQodGFnKTtcbiAgICBpZiAoY2hpbGRyZW4pIHtcbiAgICAgICAgaWYgKCEgKGNoaWxkcmVuIGluc3RhbmNlb2YgQXJyYXkpKSB7XG4gICAgICAgICAgICBjaGlsZHJlbiA9IFtjaGlsZHJlbl07XG4gICAgICAgIH1cbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBjaGlsZHJlbi5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgdmFyIGMgPSBjaGlsZHJlbltpXTtcbiAgICAgICAgICAgIGlmIChjKSB7XG4gICAgICAgICAgICAgICAgaWYgKHR5cGVvZiBjID09ICdzdHJpbmcnKSB7XG4gICAgICAgICAgICAgICAgICAgIGMgPSBkb2N1bWVudC5jcmVhdGVUZXh0Tm9kZShjKTtcbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHR5cGVvZiBjID09ICdudW1iZXInKSB7XG4gICAgICAgICAgICAgICAgICAgIGMgPSBkb2N1bWVudC5jcmVhdGVUZXh0Tm9kZSgnJyArIGMpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBlbGUuYXBwZW5kQ2hpbGQoYyk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG4gICAgXG4gICAgaWYgKGF0dHJpYnMpIHtcbiAgICAgICAgZm9yICh2YXIgbCBpbiBhdHRyaWJzKSB7XG4gICAgICAgICAgICB0cnkge1xuICAgICAgICAgICAgICAgIGVsZVtsXSA9IGF0dHJpYnNbbF07XG4gICAgICAgICAgICB9IGNhdGNoIChlKSB7XG4gICAgICAgICAgICAgICAgY29uc29sZS5sb2coJ2Vycm9yIHNldHRpbmcgJyArIGwpO1xuICAgICAgICAgICAgICAgIHRocm93KGUpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfVxuICAgIGlmIChzdHlsZXMpIHtcbiAgICAgICAgZm9yICh2YXIgbCBpbiBzdHlsZXMpIHtcbiAgICAgICAgICAgIGVsZS5zdHlsZVtsXSA9IHN0eWxlc1tsXTtcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gZWxlO1xufVxuXG5mdW5jdGlvbiBtYWtlRWxlbWVudE5TKG5hbWVzcGFjZSwgdGFnLCBjaGlsZHJlbiwgYXR0cmlicylcbntcbiAgICB2YXIgZWxlID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudE5TKG5hbWVzcGFjZSwgdGFnKTtcbiAgICBpZiAoY2hpbGRyZW4pIHtcbiAgICAgICAgaWYgKCEgKGNoaWxkcmVuIGluc3RhbmNlb2YgQXJyYXkpKSB7XG4gICAgICAgICAgICBjaGlsZHJlbiA9IFtjaGlsZHJlbl07XG4gICAgICAgIH1cbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBjaGlsZHJlbi5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgdmFyIGMgPSBjaGlsZHJlbltpXTtcbiAgICAgICAgICAgIGlmICh0eXBlb2YgYyA9PSAnc3RyaW5nJykge1xuICAgICAgICAgICAgICAgIGMgPSBkb2N1bWVudC5jcmVhdGVUZXh0Tm9kZShjKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGVsZS5hcHBlbmRDaGlsZChjKTtcbiAgICAgICAgfVxuICAgIH1cbiAgICBcbiAgICBzZXRBdHRycyhlbGUsIGF0dHJpYnMpO1xuICAgIHJldHVybiBlbGU7XG59XG5cbnZhciBhdHRyX25hbWVfY2FjaGUgPSB7fTtcblxuZnVuY3Rpb24gc2V0QXR0cihub2RlLCBrZXksIHZhbHVlKVxue1xuICAgIHZhciBhdHRyID0gYXR0cl9uYW1lX2NhY2hlW2tleV07XG4gICAgaWYgKCFhdHRyKSB7XG4gICAgICAgIHZhciBfYXR0ciA9ICcnO1xuICAgICAgICBmb3IgKHZhciBjID0gMDsgYyA8IGtleS5sZW5ndGg7ICsrYykge1xuICAgICAgICAgICAgdmFyIGNjID0ga2V5LnN1YnN0cmluZyhjLCBjKzEpO1xuICAgICAgICAgICAgdmFyIGxjYyA9IGNjLnRvTG93ZXJDYXNlKCk7XG4gICAgICAgICAgICBpZiAobGNjICE9IGNjKSB7XG4gICAgICAgICAgICAgICAgX2F0dHIgPSBfYXR0ciArICctJyArIGxjYztcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgX2F0dHIgPSBfYXR0ciArIGNjO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGF0dHJfbmFtZV9jYWNoZVtrZXldID0gX2F0dHI7XG4gICAgICAgIGF0dHIgPSBfYXR0cjtcbiAgICB9XG4gICAgbm9kZS5zZXRBdHRyaWJ1dGUoYXR0ciwgdmFsdWUpO1xufVxuXG5mdW5jdGlvbiBzZXRBdHRycyhub2RlLCBhdHRyaWJzKVxue1xuICAgIGlmIChhdHRyaWJzKSB7XG4gICAgICAgIGZvciAodmFyIGwgaW4gYXR0cmlicykge1xuICAgICAgICAgICAgc2V0QXR0cihub2RlLCBsLCBhdHRyaWJzW2xdKTtcbiAgICAgICAgfVxuICAgIH1cbn1cblxuXG5cbmZ1bmN0aW9uIHJlbW92ZUNoaWxkcmVuKG5vZGUpXG57XG4gICAgaWYgKCFub2RlIHx8ICFub2RlLmNoaWxkTm9kZXMpIHtcbiAgICAgICAgcmV0dXJuO1xuICAgIH1cblxuICAgIHdoaWxlIChub2RlLmNoaWxkTm9kZXMubGVuZ3RoID4gMCkge1xuICAgICAgICBub2RlLnJlbW92ZUNoaWxkKG5vZGUuZmlyc3RDaGlsZCk7XG4gICAgfVxufVxuXG5cblxuLy9cbi8vIFdBUk5JTkc6IG5vdCBmb3IgZ2VuZXJhbCB1c2UhXG4vL1xuXG5mdW5jdGlvbiBtaW5pSlNPTmlmeShvLCBleGMpIHtcbiAgICBpZiAodHlwZW9mIG8gPT09ICd1bmRlZmluZWQnKSB7XG4gICAgICAgIHJldHVybiAndW5kZWZpbmVkJztcbiAgICB9IGVsc2UgaWYgKG8gPT0gbnVsbCkge1xuICAgICAgICByZXR1cm4gJ251bGwnO1xuICAgIH0gZWxzZSBpZiAodHlwZW9mIG8gPT0gJ3N0cmluZycpIHtcbiAgICAgICAgcmV0dXJuIFwiJ1wiICsgbyArIFwiJ1wiO1xuICAgIH0gZWxzZSBpZiAodHlwZW9mIG8gPT0gJ251bWJlcicpIHtcbiAgICAgICAgcmV0dXJuIFwiXCIgKyBvO1xuICAgIH0gZWxzZSBpZiAodHlwZW9mIG8gPT0gJ2Jvb2xlYW4nKSB7XG4gICAgICAgIHJldHVybiBcIlwiICsgbztcbiAgICB9IGVsc2UgaWYgKHR5cGVvZiBvID09ICdvYmplY3QnKSB7XG4gICAgICAgIGlmIChvIGluc3RhbmNlb2YgQXJyYXkpIHtcbiAgICAgICAgICAgIHZhciBzID0gbnVsbDtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgby5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgICAgIHMgPSAocyA9PSBudWxsID8gJycgOiAocyArICcsICcpKSArIG1pbmlKU09OaWZ5KG9baV0sIGV4Yyk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4gJ1snICsgKHM/czonJykgKyAnXSc7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBleGMgPSBleGMgfHwge307XG4gICAgICAgICAgICB2YXIgcyA9IG51bGw7XG4gICAgICAgICAgICBmb3IgKHZhciBrIGluIG8pIHtcbiAgICAgICAgICAgICAgICBpZiAoZXhjW2tdKVxuICAgICAgICAgICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgICAgICAgICBpZiAoayAhPSB1bmRlZmluZWQgJiYgdHlwZW9mKG9ba10pICE9ICdmdW5jdGlvbicpIHtcbiAgICAgICAgICAgICAgICAgICAgcyA9IChzID09IG51bGwgPyAnJyA6IChzICsgJywgJykpICsgayArICc6ICcgKyBtaW5pSlNPTmlmeShvW2tdLCBleGMpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJldHVybiAneycgKyAocz9zOicnKSArICd9JztcbiAgICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICAgIHJldHVybiAodHlwZW9mIG8pO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gc2hhbGxvd0NvcHkobykge1xuICAgIHZhciBuID0ge307XG4gICAgZm9yICh2YXIgayBpbiBvKSB7XG4gICAgICAgIG5ba10gPSBvW2tdO1xuICAgIH1cbiAgICByZXR1cm4gbjtcbn1cblxuZnVuY3Rpb24gT2JzZXJ2ZWQoeCkge1xuICAgIHRoaXMudmFsdWUgPSB4O1xuICAgIHRoaXMubGlzdGVuZXJzID0gW107XG59XG5cbk9ic2VydmVkLnByb3RvdHlwZS5hZGRMaXN0ZW5lciA9IGZ1bmN0aW9uKGYpIHtcbiAgICB0aGlzLmxpc3RlbmVycy5wdXNoKGYpO1xufVxuXG5PYnNlcnZlZC5wcm90b3R5cGUuYWRkTGlzdGVuZXJBbmRGaXJlID0gZnVuY3Rpb24oZikge1xuICAgIHRoaXMubGlzdGVuZXJzLnB1c2goZik7XG4gICAgZih0aGlzLnZhbHVlKTtcbn1cblxuT2JzZXJ2ZWQucHJvdG90eXBlLnJlbW92ZUxpc3RlbmVyID0gZnVuY3Rpb24oZikge1xuICAgIGFycmF5UmVtb3ZlKHRoaXMubGlzdGVuZXJzLCBmKTtcbn1cblxuT2JzZXJ2ZWQucHJvdG90eXBlLmdldCA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnZhbHVlO1xufVxuXG5PYnNlcnZlZC5wcm90b3R5cGUuc2V0ID0gZnVuY3Rpb24oeCkge1xuICAgIHRoaXMudmFsdWUgPSB4O1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5saXN0ZW5lcnMubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgdGhpcy5saXN0ZW5lcnNbaV0oeCk7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBBd2FpdGVkKCkge1xuICAgIHRoaXMucXVldWUgPSBbXTtcbn1cblxuQXdhaXRlZC5wcm90b3R5cGUucHJvdmlkZSA9IGZ1bmN0aW9uKHgpIHtcbiAgICBpZiAodGhpcy5yZXMgIT09IHVuZGVmaW5lZCkge1xuICAgICAgICB0aHJvdyBcIlJlc291cmNlIGhhcyBhbHJlYWR5IGJlZW4gcHJvdmlkZWQuXCI7XG4gICAgfVxuXG4gICAgdGhpcy5yZXMgPSB4O1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5xdWV1ZS5sZW5ndGg7ICsraSkge1xuICAgICAgICB0aGlzLnF1ZXVlW2ldKHgpO1xuICAgIH1cbiAgICB0aGlzLnF1ZXVlID0gbnVsbDsgICAvLyBhdm9pZCBsZWFraW5nIGNsb3N1cmVzLlxufVxuXG5Bd2FpdGVkLnByb3RvdHlwZS5hd2FpdCA9IGZ1bmN0aW9uKGYpIHtcbiAgICBpZiAodGhpcy5yZXMgIT09IHVuZGVmaW5lZCkge1xuICAgICAgICBmKHRoaXMucmVzKTtcbiAgICAgICAgcmV0dXJuIHRoaXMucmVzO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHRoaXMucXVldWUucHVzaChmKTtcbiAgICB9XG59XG5cbnZhciBfX2RhbGxpYW5jZV9zYWx0U2VlZCA9IDA7XG5cbmZ1bmN0aW9uIHNhbHRVUkwodXJsKSB7XG4gICAgcmV0dXJuIHVybCArICc/c2FsdD0nICsgYjY0X3NoYTEoJycgKyBEYXRlLm5vdygpICsgJywnICsgKCsrX19kYWxsaWFuY2Vfc2FsdFNlZWQpKTtcbn1cblxuZnVuY3Rpb24gdGV4dFhIUih1cmwsIGNhbGxiYWNrLCBvcHRzKSB7XG4gICAgaWYgKG9wdHMuc2FsdCkgXG4gICAgICAgIHVybCA9IHNhbHRVUkwodXJsKTtcblxuICAgIHZhciByZXEgPSBuZXcgWE1MSHR0cFJlcXVlc3QoKTtcbiAgICByZXEub25yZWFkeXN0YXRlY2hhbmdlID0gZnVuY3Rpb24oKSB7XG4gICAgXHRpZiAocmVxLnJlYWR5U3RhdGUgPT0gNCkge1xuICAgIFx0ICAgIGlmIChyZXEuc3RhdHVzID49IDMwMCkge1xuICAgIFx0XHQgICAgY2FsbGJhY2sobnVsbCwgJ0Vycm9yIGNvZGUgJyArIHJlcS5zdGF0dXMpO1xuICAgIFx0ICAgIH0gZWxzZSB7XG4gICAgXHRcdCAgICBjYWxsYmFjayhyZXEucmVzcG9uc2VUZXh0KTtcbiAgICBcdCAgICB9XG4gICAgXHR9XG4gICAgfTtcbiAgICBcbiAgICByZXEub3BlbignR0VUJywgdXJsLCB0cnVlKTtcbiAgICByZXEucmVzcG9uc2VUeXBlID0gJ3RleHQnO1xuXG4gICAgaWYgKG9wdHMgJiYgb3B0cy5jcmVkZW50aWFscykge1xuICAgICAgICByZXEud2l0aENyZWRlbnRpYWxzID0gdHJ1ZTtcbiAgICB9XG4gICAgcmVxLnNlbmQoJycpO1xufVxuXG5mdW5jdGlvbiByZWxhdGl2ZVVSTChiYXNlLCByZWwpIHtcbiAgICAvLyBGSVhNRSBxdWl0ZSBuYWl2ZSAtLSBnb29kIGVub3VnaCBmb3IgdHJhY2todWJzP1xuXG4gICAgaWYgKHJlbC5pbmRleE9mKCdodHRwOicpID09IDAgfHwgcmVsLmluZGV4T2YoJ2h0dHBzOicpID09IDApIHtcbiAgICAgICAgcmV0dXJuIHJlbDtcbiAgICB9XG5cbiAgICB2YXIgbGkgPSBiYXNlLmxhc3RJbmRleE9mKCcvJyk7XG4gICAgaWYgKGxpID49IDApIHtcbiAgICAgICAgcmV0dXJuIGJhc2Uuc3Vic3RyKDAsIGxpICsgMSkgKyByZWw7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgcmV0dXJuIHJlbDtcbiAgICB9XG59XG5cbnZhciBBTUlOT19BQ0lEX1RSQU5TTEFUSU9OID0ge1xuICAgICdUVFQnOiAnRicsXG4gICAgJ1RUQyc6ICdGJyxcbiAgICAnVFRBJzogJ0wnLFxuICAgICdUVEcnOiAnTCcsXG4gICAgJ0NUVCc6ICdMJyxcbiAgICAnQ1RDJzogJ0wnLFxuICAgICdDVEEnOiAnTCcsXG4gICAgJ0NURyc6ICdMJyxcbiAgICAnQVRUJzogJ0knLFxuICAgICdBVEMnOiAnSScsXG4gICAgJ0FUQSc6ICdJJyxcbiAgICAnQVRHJzogJ00nLFxuICAgICdHVFQnOiAnVicsXG4gICAgJ0dUQyc6ICdWJyxcbiAgICAnR1RBJzogJ1YnLFxuICAgICdHVEcnOiAnVicsXG4gICAgJ1RDVCc6ICdTJyxcbiAgICAnVENDJzogJ1MnLFxuICAgICdUQ0EnOiAnUycsXG4gICAgJ1RDRyc6ICdTJyxcbiAgICAnQ0NUJzogJ1AnLFxuICAgICdDQ0MnOiAnUCcsXG4gICAgJ0NDQSc6ICdQJyxcbiAgICAnQ0NHJzogJ1AnLFxuICAgICdBQ1QnOiAnVCcsXG4gICAgJ0FDQyc6ICdUJyxcbiAgICAnQUNBJzogJ1QnLFxuICAgICdBQ0cnOiAnVCcsXG4gICAgJ0dDVCc6ICdBJyxcbiAgICAnR0NDJzogJ0EnLFxuICAgICdHQ0EnOiAnQScsXG4gICAgJ0dDRyc6ICdBJyxcbiAgICAnVEFUJzogJ1knLFxuICAgICdUQUMnOiAnWScsXG4gICAgJ1RBQSc6ICcqJywgIC8vIHN0b3BcbiAgICAnVEFHJzogJyonLCAgLy8gc3RvcFxuICAgICdDQVQnOiAnSCcsXG4gICAgJ0NBQyc6ICdIJyxcbiAgICAnQ0FBJzogJ1EnLFxuICAgICdDQUcnOiAnUScsXG4gICAgJ0FBVCc6ICdOJyxcbiAgICAnQUFDJzogJ04nLFxuICAgICdBQUEnOiAnSycsXG4gICAgJ0FBRyc6ICdLJyxcbiAgICAnR0FUJzogJ0QnLFxuICAgICdHQUMnOiAnRCcsXG4gICAgJ0dBQSc6ICdFJyxcbiAgICAnR0FHJzogJ0UnLFxuICAgICdUR1QnOiAnQycsXG4gICAgJ1RHQyc6ICdDJyxcbiAgICAnVEdBJzogJyonLCAgLy8gc3RvcFxuICAgICdUR0cnOiAnVycsXG4gICAgJ0NHVCc6ICdSJyxcbiAgICAnQ0dDJzogJ1InLFxuICAgICdDR0EnOiAnUicsXG4gICAgJ0NHRyc6ICdSJyxcbiAgICAnQUdUJzogJ1MnLFxuICAgICdBR0MnOiAnUycsXG4gICAgJ0FHQSc6ICdSJyxcbiAgICAnQUdHJzogJ1InLFxuICAgICdHR1QnOiAnRycsXG4gICAgJ0dHQyc6ICdHJyxcbiAgICAnR0dBJzogJ0cnLFxuICAgICdHR0cnOiAnRydcbn1cblxuZnVuY3Rpb24gcmVzb2x2ZVVybFRvUGFnZShyZWwpIHtcbiAgICByZXR1cm4gbWFrZUVsZW1lbnQoJ2EnLCBudWxsLCB7aHJlZjogcmVsfSkuaHJlZjtcbn1cblxuLy9cbi8vIE1pc3NpbmcgQVBJc1xuLy8gXG5cbmlmICghKCd0cmltJyBpbiBTdHJpbmcucHJvdG90eXBlKSkge1xuICAgIFN0cmluZy5wcm90b3R5cGUudHJpbSA9IGZ1bmN0aW9uKCkge1xuICAgICAgICByZXR1cm4gdGhpcy5yZXBsYWNlKC9eXFxzKy8sICcnKS5yZXBsYWNlKC9cXHMrJC8sICcnKTtcbiAgICB9O1xufVxuXG5pZiAodHlwZW9mKG1vZHVsZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgICAgIHRleHRYSFI6IHRleHRYSFIsXG4gICAgICAgIHJlbGF0aXZlVVJMOiByZWxhdGl2ZVVSTCxcbiAgICAgICAgcmVzb2x2ZVVybFRvUGFnZTogcmVzb2x2ZVVybFRvUGFnZSxcbiAgICAgICAgc2hhbGxvd0NvcHk6IHNoYWxsb3dDb3B5LFxuICAgICAgICBwdXNobzogcHVzaG8sXG4gICAgICAgIHB1c2huZXc6IHB1c2huZXcsXG4gICAgICAgIHB1c2huZXdvOiBwdXNobmV3byxcbiAgICAgICAgYXJyYXlJbmRleE9mOiBhcnJheUluZGV4T2YsXG4gICAgICAgIHBpY2s6IHBpY2ssXG5cbiAgICAgICAgbWFrZUVsZW1lbnQ6IG1ha2VFbGVtZW50LFxuICAgICAgICBtYWtlRWxlbWVudE5TOiBtYWtlRWxlbWVudE5TLFxuICAgICAgICByZW1vdmVDaGlsZHJlbjogcmVtb3ZlQ2hpbGRyZW4sXG5cbiAgICAgICAgbWluaUpTT05pZnk6IG1pbmlKU09OaWZ5LFxuXG4gICAgICAgIE9ic2VydmVkOiBPYnNlcnZlZCxcbiAgICAgICAgQXdhaXRlZDogQXdhaXRlZCxcblxuICAgICAgICBBTUlOT19BQ0lEX1RSQU5TTEFUSU9OOiBBTUlOT19BQ0lEX1RSQU5TTEFUSU9OXG4gICAgfVxufVxuIiwiLyogLSotIG1vZGU6IGphdmFzY3JpcHQ7IGMtYmFzaWMtb2Zmc2V0OiA0OyBpbmRlbnQtdGFicy1tb2RlOiBuaWwgLSotICovXG5cbi8vIFxuLy8gSmF2YXNjcmlwdCBaTGliXG4vLyBCeSBUaG9tYXMgRG93biAyMDEwLTIwMTFcbi8vXG4vLyBCYXNlZCB2ZXJ5IGhlYXZpbHkgb24gcG9ydGlvbnMgb2YganpsaWIgKGJ5IHltbmtAamNyYWZ0LmNvbSksIHdobyBpblxuLy8gdHVybiBjcmVkaXRzIEplYW4tbG91cCBHYWlsbHkgYW5kIE1hcmsgQWRsZXIgZm9yIHRoZSBvcmlnaW5hbCB6bGliIGNvZGUuXG4vL1xuLy8gaW5mbGF0ZS5qczogWkxpYiBpbmZsYXRlIGNvZGVcbi8vXG5cbi8vXG4vLyBTaGFyZWQgY29uc3RhbnRzXG4vL1xuXG52YXIgTUFYX1dCSVRTPTE1OyAvLyAzMksgTFo3NyB3aW5kb3dcbnZhciBERUZfV0JJVFM9TUFYX1dCSVRTO1xudmFyIE1BWF9NRU1fTEVWRUw9OTtcbnZhciBNQU5ZPTE0NDA7XG52YXIgQk1BWCA9IDE1O1xuXG4vLyBwcmVzZXQgZGljdGlvbmFyeSBmbGFnIGluIHpsaWIgaGVhZGVyXG52YXIgUFJFU0VUX0RJQ1Q9MHgyMDtcblxudmFyIFpfTk9fRkxVU0g9MDtcbnZhciBaX1BBUlRJQUxfRkxVU0g9MTtcbnZhciBaX1NZTkNfRkxVU0g9MjtcbnZhciBaX0ZVTExfRkxVU0g9MztcbnZhciBaX0ZJTklTSD00O1xuXG52YXIgWl9ERUZMQVRFRD04O1xuXG52YXIgWl9PSz0wO1xudmFyIFpfU1RSRUFNX0VORD0xO1xudmFyIFpfTkVFRF9ESUNUPTI7XG52YXIgWl9FUlJOTz0tMTtcbnZhciBaX1NUUkVBTV9FUlJPUj0tMjtcbnZhciBaX0RBVEFfRVJST1I9LTM7XG52YXIgWl9NRU1fRVJST1I9LTQ7XG52YXIgWl9CVUZfRVJST1I9LTU7XG52YXIgWl9WRVJTSU9OX0VSUk9SPS02O1xuXG52YXIgTUVUSE9EPTA7ICAgLy8gd2FpdGluZyBmb3IgbWV0aG9kIGJ5dGVcbnZhciBGTEFHPTE7ICAgICAvLyB3YWl0aW5nIGZvciBmbGFnIGJ5dGVcbnZhciBESUNUND0yOyAgICAvLyBmb3VyIGRpY3Rpb25hcnkgY2hlY2sgYnl0ZXMgdG8gZ29cbnZhciBESUNUMz0zOyAgICAvLyB0aHJlZSBkaWN0aW9uYXJ5IGNoZWNrIGJ5dGVzIHRvIGdvXG52YXIgRElDVDI9NDsgICAgLy8gdHdvIGRpY3Rpb25hcnkgY2hlY2sgYnl0ZXMgdG8gZ29cbnZhciBESUNUMT01OyAgICAvLyBvbmUgZGljdGlvbmFyeSBjaGVjayBieXRlIHRvIGdvXG52YXIgRElDVDA9NjsgICAgLy8gd2FpdGluZyBmb3IgaW5mbGF0ZVNldERpY3Rpb25hcnlcbnZhciBCTE9DS1M9NzsgICAvLyBkZWNvbXByZXNzaW5nIGJsb2Nrc1xudmFyIENIRUNLND04OyAgIC8vIGZvdXIgY2hlY2sgYnl0ZXMgdG8gZ29cbnZhciBDSEVDSzM9OTsgICAvLyB0aHJlZSBjaGVjayBieXRlcyB0byBnb1xudmFyIENIRUNLMj0xMDsgIC8vIHR3byBjaGVjayBieXRlcyB0byBnb1xudmFyIENIRUNLMT0xMTsgIC8vIG9uZSBjaGVjayBieXRlIHRvIGdvXG52YXIgRE9ORT0xMjsgICAgLy8gZmluaXNoZWQgY2hlY2ssIGRvbmVcbnZhciBCQUQ9MTM7ICAgICAvLyBnb3QgYW4gZXJyb3ItLXN0YXkgaGVyZVxuXG52YXIgaW5mbGF0ZV9tYXNrID0gWzB4MDAwMDAwMDAsIDB4MDAwMDAwMDEsIDB4MDAwMDAwMDMsIDB4MDAwMDAwMDcsIDB4MDAwMDAwMGYsIDB4MDAwMDAwMWYsIDB4MDAwMDAwM2YsIDB4MDAwMDAwN2YsIDB4MDAwMDAwZmYsIDB4MDAwMDAxZmYsIDB4MDAwMDAzZmYsIDB4MDAwMDA3ZmYsIDB4MDAwMDBmZmYsIDB4MDAwMDFmZmYsIDB4MDAwMDNmZmYsIDB4MDAwMDdmZmYsIDB4MDAwMGZmZmZdO1xuXG52YXIgSUJfVFlQRT0wOyAgLy8gZ2V0IHR5cGUgYml0cyAoMywgaW5jbHVkaW5nIGVuZCBiaXQpXG52YXIgSUJfTEVOUz0xOyAgLy8gZ2V0IGxlbmd0aHMgZm9yIHN0b3JlZFxudmFyIElCX1NUT1JFRD0yOy8vIHByb2Nlc3Npbmcgc3RvcmVkIGJsb2NrXG52YXIgSUJfVEFCTEU9MzsgLy8gZ2V0IHRhYmxlIGxlbmd0aHNcbnZhciBJQl9CVFJFRT00OyAvLyBnZXQgYml0IGxlbmd0aHMgdHJlZSBmb3IgYSBkeW5hbWljIGJsb2NrXG52YXIgSUJfRFRSRUU9NTsgLy8gZ2V0IGxlbmd0aCwgZGlzdGFuY2UgdHJlZXMgZm9yIGEgZHluYW1pYyBibG9ja1xudmFyIElCX0NPREVTPTY7IC8vIHByb2Nlc3NpbmcgZml4ZWQgb3IgZHluYW1pYyBibG9ja1xudmFyIElCX0RSWT03OyAgIC8vIG91dHB1dCByZW1haW5pbmcgd2luZG93IGJ5dGVzXG52YXIgSUJfRE9ORT04OyAgLy8gZmluaXNoZWQgbGFzdCBibG9jaywgZG9uZVxudmFyIElCX0JBRD05OyAgIC8vIG90IGEgZGF0YSBlcnJvci0tc3R1Y2sgaGVyZVxuXG52YXIgZml4ZWRfYmwgPSA5O1xudmFyIGZpeGVkX2JkID0gNTtcblxudmFyIGZpeGVkX3RsID0gW1xuICAgIDk2LDcsMjU2LCAwLDgsODAsIDAsOCwxNiwgODQsOCwxMTUsXG4gICAgODIsNywzMSwgMCw4LDExMiwgMCw4LDQ4LCAwLDksMTkyLFxuICAgIDgwLDcsMTAsIDAsOCw5NiwgMCw4LDMyLCAwLDksMTYwLFxuICAgIDAsOCwwLCAwLDgsMTI4LCAwLDgsNjQsIDAsOSwyMjQsXG4gICAgODAsNyw2LCAwLDgsODgsIDAsOCwyNCwgMCw5LDE0NCxcbiAgICA4Myw3LDU5LCAwLDgsMTIwLCAwLDgsNTYsIDAsOSwyMDgsXG4gICAgODEsNywxNywgMCw4LDEwNCwgMCw4LDQwLCAwLDksMTc2LFxuICAgIDAsOCw4LCAwLDgsMTM2LCAwLDgsNzIsIDAsOSwyNDAsXG4gICAgODAsNyw0LCAwLDgsODQsIDAsOCwyMCwgODUsOCwyMjcsXG4gICAgODMsNyw0MywgMCw4LDExNiwgMCw4LDUyLCAwLDksMjAwLFxuICAgIDgxLDcsMTMsIDAsOCwxMDAsIDAsOCwzNiwgMCw5LDE2OCxcbiAgICAwLDgsNCwgMCw4LDEzMiwgMCw4LDY4LCAwLDksMjMyLFxuICAgIDgwLDcsOCwgMCw4LDkyLCAwLDgsMjgsIDAsOSwxNTIsXG4gICAgODQsNyw4MywgMCw4LDEyNCwgMCw4LDYwLCAwLDksMjE2LFxuICAgIDgyLDcsMjMsIDAsOCwxMDgsIDAsOCw0NCwgMCw5LDE4NCxcbiAgICAwLDgsMTIsIDAsOCwxNDAsIDAsOCw3NiwgMCw5LDI0OCxcbiAgICA4MCw3LDMsIDAsOCw4MiwgMCw4LDE4LCA4NSw4LDE2MyxcbiAgICA4Myw3LDM1LCAwLDgsMTE0LCAwLDgsNTAsIDAsOSwxOTYsXG4gICAgODEsNywxMSwgMCw4LDk4LCAwLDgsMzQsIDAsOSwxNjQsXG4gICAgMCw4LDIsIDAsOCwxMzAsIDAsOCw2NiwgMCw5LDIyOCxcbiAgICA4MCw3LDcsIDAsOCw5MCwgMCw4LDI2LCAwLDksMTQ4LFxuICAgIDg0LDcsNjcsIDAsOCwxMjIsIDAsOCw1OCwgMCw5LDIxMixcbiAgICA4Miw3LDE5LCAwLDgsMTA2LCAwLDgsNDIsIDAsOSwxODAsXG4gICAgMCw4LDEwLCAwLDgsMTM4LCAwLDgsNzQsIDAsOSwyNDQsXG4gICAgODAsNyw1LCAwLDgsODYsIDAsOCwyMiwgMTkyLDgsMCxcbiAgICA4Myw3LDUxLCAwLDgsMTE4LCAwLDgsNTQsIDAsOSwyMDQsXG4gICAgODEsNywxNSwgMCw4LDEwMiwgMCw4LDM4LCAwLDksMTcyLFxuICAgIDAsOCw2LCAwLDgsMTM0LCAwLDgsNzAsIDAsOSwyMzYsXG4gICAgODAsNyw5LCAwLDgsOTQsIDAsOCwzMCwgMCw5LDE1NixcbiAgICA4NCw3LDk5LCAwLDgsMTI2LCAwLDgsNjIsIDAsOSwyMjAsXG4gICAgODIsNywyNywgMCw4LDExMCwgMCw4LDQ2LCAwLDksMTg4LFxuICAgIDAsOCwxNCwgMCw4LDE0MiwgMCw4LDc4LCAwLDksMjUyLFxuICAgIDk2LDcsMjU2LCAwLDgsODEsIDAsOCwxNywgODUsOCwxMzEsXG4gICAgODIsNywzMSwgMCw4LDExMywgMCw4LDQ5LCAwLDksMTk0LFxuICAgIDgwLDcsMTAsIDAsOCw5NywgMCw4LDMzLCAwLDksMTYyLFxuICAgIDAsOCwxLCAwLDgsMTI5LCAwLDgsNjUsIDAsOSwyMjYsXG4gICAgODAsNyw2LCAwLDgsODksIDAsOCwyNSwgMCw5LDE0NixcbiAgICA4Myw3LDU5LCAwLDgsMTIxLCAwLDgsNTcsIDAsOSwyMTAsXG4gICAgODEsNywxNywgMCw4LDEwNSwgMCw4LDQxLCAwLDksMTc4LFxuICAgIDAsOCw5LCAwLDgsMTM3LCAwLDgsNzMsIDAsOSwyNDIsXG4gICAgODAsNyw0LCAwLDgsODUsIDAsOCwyMSwgODAsOCwyNTgsXG4gICAgODMsNyw0MywgMCw4LDExNywgMCw4LDUzLCAwLDksMjAyLFxuICAgIDgxLDcsMTMsIDAsOCwxMDEsIDAsOCwzNywgMCw5LDE3MCxcbiAgICAwLDgsNSwgMCw4LDEzMywgMCw4LDY5LCAwLDksMjM0LFxuICAgIDgwLDcsOCwgMCw4LDkzLCAwLDgsMjksIDAsOSwxNTQsXG4gICAgODQsNyw4MywgMCw4LDEyNSwgMCw4LDYxLCAwLDksMjE4LFxuICAgIDgyLDcsMjMsIDAsOCwxMDksIDAsOCw0NSwgMCw5LDE4NixcbiAgICAwLDgsMTMsIDAsOCwxNDEsIDAsOCw3NywgMCw5LDI1MCxcbiAgICA4MCw3LDMsIDAsOCw4MywgMCw4LDE5LCA4NSw4LDE5NSxcbiAgICA4Myw3LDM1LCAwLDgsMTE1LCAwLDgsNTEsIDAsOSwxOTgsXG4gICAgODEsNywxMSwgMCw4LDk5LCAwLDgsMzUsIDAsOSwxNjYsXG4gICAgMCw4LDMsIDAsOCwxMzEsIDAsOCw2NywgMCw5LDIzMCxcbiAgICA4MCw3LDcsIDAsOCw5MSwgMCw4LDI3LCAwLDksMTUwLFxuICAgIDg0LDcsNjcsIDAsOCwxMjMsIDAsOCw1OSwgMCw5LDIxNCxcbiAgICA4Miw3LDE5LCAwLDgsMTA3LCAwLDgsNDMsIDAsOSwxODIsXG4gICAgMCw4LDExLCAwLDgsMTM5LCAwLDgsNzUsIDAsOSwyNDYsXG4gICAgODAsNyw1LCAwLDgsODcsIDAsOCwyMywgMTkyLDgsMCxcbiAgICA4Myw3LDUxLCAwLDgsMTE5LCAwLDgsNTUsIDAsOSwyMDYsXG4gICAgODEsNywxNSwgMCw4LDEwMywgMCw4LDM5LCAwLDksMTc0LFxuICAgIDAsOCw3LCAwLDgsMTM1LCAwLDgsNzEsIDAsOSwyMzgsXG4gICAgODAsNyw5LCAwLDgsOTUsIDAsOCwzMSwgMCw5LDE1OCxcbiAgICA4NCw3LDk5LCAwLDgsMTI3LCAwLDgsNjMsIDAsOSwyMjIsXG4gICAgODIsNywyNywgMCw4LDExMSwgMCw4LDQ3LCAwLDksMTkwLFxuICAgIDAsOCwxNSwgMCw4LDE0MywgMCw4LDc5LCAwLDksMjU0LFxuICAgIDk2LDcsMjU2LCAwLDgsODAsIDAsOCwxNiwgODQsOCwxMTUsXG4gICAgODIsNywzMSwgMCw4LDExMiwgMCw4LDQ4LCAwLDksMTkzLFxuXG4gICAgODAsNywxMCwgMCw4LDk2LCAwLDgsMzIsIDAsOSwxNjEsXG4gICAgMCw4LDAsIDAsOCwxMjgsIDAsOCw2NCwgMCw5LDIyNSxcbiAgICA4MCw3LDYsIDAsOCw4OCwgMCw4LDI0LCAwLDksMTQ1LFxuICAgIDgzLDcsNTksIDAsOCwxMjAsIDAsOCw1NiwgMCw5LDIwOSxcbiAgICA4MSw3LDE3LCAwLDgsMTA0LCAwLDgsNDAsIDAsOSwxNzcsXG4gICAgMCw4LDgsIDAsOCwxMzYsIDAsOCw3MiwgMCw5LDI0MSxcbiAgICA4MCw3LDQsIDAsOCw4NCwgMCw4LDIwLCA4NSw4LDIyNyxcbiAgICA4Myw3LDQzLCAwLDgsMTE2LCAwLDgsNTIsIDAsOSwyMDEsXG4gICAgODEsNywxMywgMCw4LDEwMCwgMCw4LDM2LCAwLDksMTY5LFxuICAgIDAsOCw0LCAwLDgsMTMyLCAwLDgsNjgsIDAsOSwyMzMsXG4gICAgODAsNyw4LCAwLDgsOTIsIDAsOCwyOCwgMCw5LDE1MyxcbiAgICA4NCw3LDgzLCAwLDgsMTI0LCAwLDgsNjAsIDAsOSwyMTcsXG4gICAgODIsNywyMywgMCw4LDEwOCwgMCw4LDQ0LCAwLDksMTg1LFxuICAgIDAsOCwxMiwgMCw4LDE0MCwgMCw4LDc2LCAwLDksMjQ5LFxuICAgIDgwLDcsMywgMCw4LDgyLCAwLDgsMTgsIDg1LDgsMTYzLFxuICAgIDgzLDcsMzUsIDAsOCwxMTQsIDAsOCw1MCwgMCw5LDE5NyxcbiAgICA4MSw3LDExLCAwLDgsOTgsIDAsOCwzNCwgMCw5LDE2NSxcbiAgICAwLDgsMiwgMCw4LDEzMCwgMCw4LDY2LCAwLDksMjI5LFxuICAgIDgwLDcsNywgMCw4LDkwLCAwLDgsMjYsIDAsOSwxNDksXG4gICAgODQsNyw2NywgMCw4LDEyMiwgMCw4LDU4LCAwLDksMjEzLFxuICAgIDgyLDcsMTksIDAsOCwxMDYsIDAsOCw0MiwgMCw5LDE4MSxcbiAgICAwLDgsMTAsIDAsOCwxMzgsIDAsOCw3NCwgMCw5LDI0NSxcbiAgICA4MCw3LDUsIDAsOCw4NiwgMCw4LDIyLCAxOTIsOCwwLFxuICAgIDgzLDcsNTEsIDAsOCwxMTgsIDAsOCw1NCwgMCw5LDIwNSxcbiAgICA4MSw3LDE1LCAwLDgsMTAyLCAwLDgsMzgsIDAsOSwxNzMsXG4gICAgMCw4LDYsIDAsOCwxMzQsIDAsOCw3MCwgMCw5LDIzNyxcbiAgICA4MCw3LDksIDAsOCw5NCwgMCw4LDMwLCAwLDksMTU3LFxuICAgIDg0LDcsOTksIDAsOCwxMjYsIDAsOCw2MiwgMCw5LDIyMSxcbiAgICA4Miw3LDI3LCAwLDgsMTEwLCAwLDgsNDYsIDAsOSwxODksXG4gICAgMCw4LDE0LCAwLDgsMTQyLCAwLDgsNzgsIDAsOSwyNTMsXG4gICAgOTYsNywyNTYsIDAsOCw4MSwgMCw4LDE3LCA4NSw4LDEzMSxcbiAgICA4Miw3LDMxLCAwLDgsMTEzLCAwLDgsNDksIDAsOSwxOTUsXG4gICAgODAsNywxMCwgMCw4LDk3LCAwLDgsMzMsIDAsOSwxNjMsXG4gICAgMCw4LDEsIDAsOCwxMjksIDAsOCw2NSwgMCw5LDIyNyxcbiAgICA4MCw3LDYsIDAsOCw4OSwgMCw4LDI1LCAwLDksMTQ3LFxuICAgIDgzLDcsNTksIDAsOCwxMjEsIDAsOCw1NywgMCw5LDIxMSxcbiAgICA4MSw3LDE3LCAwLDgsMTA1LCAwLDgsNDEsIDAsOSwxNzksXG4gICAgMCw4LDksIDAsOCwxMzcsIDAsOCw3MywgMCw5LDI0MyxcbiAgICA4MCw3LDQsIDAsOCw4NSwgMCw4LDIxLCA4MCw4LDI1OCxcbiAgICA4Myw3LDQzLCAwLDgsMTE3LCAwLDgsNTMsIDAsOSwyMDMsXG4gICAgODEsNywxMywgMCw4LDEwMSwgMCw4LDM3LCAwLDksMTcxLFxuICAgIDAsOCw1LCAwLDgsMTMzLCAwLDgsNjksIDAsOSwyMzUsXG4gICAgODAsNyw4LCAwLDgsOTMsIDAsOCwyOSwgMCw5LDE1NSxcbiAgICA4NCw3LDgzLCAwLDgsMTI1LCAwLDgsNjEsIDAsOSwyMTksXG4gICAgODIsNywyMywgMCw4LDEwOSwgMCw4LDQ1LCAwLDksMTg3LFxuICAgIDAsOCwxMywgMCw4LDE0MSwgMCw4LDc3LCAwLDksMjUxLFxuICAgIDgwLDcsMywgMCw4LDgzLCAwLDgsMTksIDg1LDgsMTk1LFxuICAgIDgzLDcsMzUsIDAsOCwxMTUsIDAsOCw1MSwgMCw5LDE5OSxcbiAgICA4MSw3LDExLCAwLDgsOTksIDAsOCwzNSwgMCw5LDE2NyxcbiAgICAwLDgsMywgMCw4LDEzMSwgMCw4LDY3LCAwLDksMjMxLFxuICAgIDgwLDcsNywgMCw4LDkxLCAwLDgsMjcsIDAsOSwxNTEsXG4gICAgODQsNyw2NywgMCw4LDEyMywgMCw4LDU5LCAwLDksMjE1LFxuICAgIDgyLDcsMTksIDAsOCwxMDcsIDAsOCw0MywgMCw5LDE4MyxcbiAgICAwLDgsMTEsIDAsOCwxMzksIDAsOCw3NSwgMCw5LDI0NyxcbiAgICA4MCw3LDUsIDAsOCw4NywgMCw4LDIzLCAxOTIsOCwwLFxuICAgIDgzLDcsNTEsIDAsOCwxMTksIDAsOCw1NSwgMCw5LDIwNyxcbiAgICA4MSw3LDE1LCAwLDgsMTAzLCAwLDgsMzksIDAsOSwxNzUsXG4gICAgMCw4LDcsIDAsOCwxMzUsIDAsOCw3MSwgMCw5LDIzOSxcbiAgICA4MCw3LDksIDAsOCw5NSwgMCw4LDMxLCAwLDksMTU5LFxuICAgIDg0LDcsOTksIDAsOCwxMjcsIDAsOCw2MywgMCw5LDIyMyxcbiAgICA4Miw3LDI3LCAwLDgsMTExLCAwLDgsNDcsIDAsOSwxOTEsXG4gICAgMCw4LDE1LCAwLDgsMTQzLCAwLDgsNzksIDAsOSwyNTVcbl07XG52YXIgZml4ZWRfdGQgPSBbXG4gICAgODAsNSwxLCA4Nyw1LDI1NywgODMsNSwxNywgOTEsNSw0MDk3LFxuICAgIDgxLDUsNSwgODksNSwxMDI1LCA4NSw1LDY1LCA5Myw1LDE2Mzg1LFxuICAgIDgwLDUsMywgODgsNSw1MTMsIDg0LDUsMzMsIDkyLDUsODE5MyxcbiAgICA4Miw1LDksIDkwLDUsMjA0OSwgODYsNSwxMjksIDE5Miw1LDI0NTc3LFxuICAgIDgwLDUsMiwgODcsNSwzODUsIDgzLDUsMjUsIDkxLDUsNjE0NSxcbiAgICA4MSw1LDcsIDg5LDUsMTUzNywgODUsNSw5NywgOTMsNSwyNDU3NyxcbiAgICA4MCw1LDQsIDg4LDUsNzY5LCA4NCw1LDQ5LCA5Miw1LDEyMjg5LFxuICAgIDgyLDUsMTMsIDkwLDUsMzA3MywgODYsNSwxOTMsIDE5Miw1LDI0NTc3XG5dO1xuXG4gIC8vIFRhYmxlcyBmb3IgZGVmbGF0ZSBmcm9tIFBLWklQJ3MgYXBwbm90ZS50eHQuXG4gIHZhciBjcGxlbnMgPSBbIC8vIENvcHkgbGVuZ3RocyBmb3IgbGl0ZXJhbCBjb2RlcyAyNTcuLjI4NVxuICAgICAgICAzLCA0LCA1LCA2LCA3LCA4LCA5LCAxMCwgMTEsIDEzLCAxNSwgMTcsIDE5LCAyMywgMjcsIDMxLFxuICAgICAgICAzNSwgNDMsIDUxLCA1OSwgNjcsIDgzLCA5OSwgMTE1LCAxMzEsIDE2MywgMTk1LCAyMjcsIDI1OCwgMCwgMFxuICBdO1xuXG4gIC8vIHNlZSBub3RlICMxMyBhYm92ZSBhYm91dCAyNThcbiAgdmFyIGNwbGV4dCA9IFsgLy8gRXh0cmEgYml0cyBmb3IgbGl0ZXJhbCBjb2RlcyAyNTcuLjI4NVxuICAgICAgICAwLCAwLCAwLCAwLCAwLCAwLCAwLCAwLCAxLCAxLCAxLCAxLCAyLCAyLCAyLCAyLFxuICAgICAgICAzLCAzLCAzLCAzLCA0LCA0LCA0LCA0LCA1LCA1LCA1LCA1LCAwLCAxMTIsIDExMiAgLy8gMTEyPT1pbnZhbGlkXG4gIF07XG5cbiB2YXIgY3BkaXN0ID0gWyAvLyBDb3B5IG9mZnNldHMgZm9yIGRpc3RhbmNlIGNvZGVzIDAuLjI5XG4gICAgICAgIDEsIDIsIDMsIDQsIDUsIDcsIDksIDEzLCAxNywgMjUsIDMzLCA0OSwgNjUsIDk3LCAxMjksIDE5MyxcbiAgICAgICAgMjU3LCAzODUsIDUxMywgNzY5LCAxMDI1LCAxNTM3LCAyMDQ5LCAzMDczLCA0MDk3LCA2MTQ1LFxuICAgICAgICA4MTkzLCAxMjI4OSwgMTYzODUsIDI0NTc3XG4gIF07XG5cbiAgdmFyIGNwZGV4dCA9IFsgLy8gRXh0cmEgYml0cyBmb3IgZGlzdGFuY2UgY29kZXNcbiAgICAgICAgMCwgMCwgMCwgMCwgMSwgMSwgMiwgMiwgMywgMywgNCwgNCwgNSwgNSwgNiwgNixcbiAgICAgICAgNywgNywgOCwgOCwgOSwgOSwgMTAsIDEwLCAxMSwgMTEsXG4gICAgICAgIDEyLCAxMiwgMTMsIDEzXTtcblxuLy9cbi8vIFpTdHJlYW0uamF2YVxuLy9cblxuZnVuY3Rpb24gWlN0cmVhbSgpIHtcbn1cblxuXG5aU3RyZWFtLnByb3RvdHlwZS5pbmZsYXRlSW5pdCA9IGZ1bmN0aW9uKHcsIG5vd3JhcCkge1xuICAgIGlmICghdykge1xuXHR3ID0gREVGX1dCSVRTO1xuICAgIH1cbiAgICBpZiAobm93cmFwKSB7XG5cdG5vd3JhcCA9IGZhbHNlO1xuICAgIH1cbiAgICB0aGlzLmlzdGF0ZSA9IG5ldyBJbmZsYXRlKCk7XG4gICAgcmV0dXJuIHRoaXMuaXN0YXRlLmluZmxhdGVJbml0KHRoaXMsIG5vd3JhcD8tdzp3KTtcbn1cblxuWlN0cmVhbS5wcm90b3R5cGUuaW5mbGF0ZSA9IGZ1bmN0aW9uKGYpIHtcbiAgICBpZih0aGlzLmlzdGF0ZT09bnVsbCkgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIHJldHVybiB0aGlzLmlzdGF0ZS5pbmZsYXRlKHRoaXMsIGYpO1xufVxuXG5aU3RyZWFtLnByb3RvdHlwZS5pbmZsYXRlRW5kID0gZnVuY3Rpb24oKXtcbiAgICBpZih0aGlzLmlzdGF0ZT09bnVsbCkgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIHZhciByZXQ9aXN0YXRlLmluZmxhdGVFbmQodGhpcyk7XG4gICAgdGhpcy5pc3RhdGUgPSBudWxsO1xuICAgIHJldHVybiByZXQ7XG59XG5aU3RyZWFtLnByb3RvdHlwZS5pbmZsYXRlU3luYyA9IGZ1bmN0aW9uKCl7XG4gICAgLy8gaWYoaXN0YXRlID09IG51bGwpIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICByZXR1cm4gaXN0YXRlLmluZmxhdGVTeW5jKHRoaXMpO1xufVxuWlN0cmVhbS5wcm90b3R5cGUuaW5mbGF0ZVNldERpY3Rpb25hcnkgPSBmdW5jdGlvbihkaWN0aW9uYXJ5LCBkaWN0TGVuZ3RoKXtcbiAgICAvLyBpZihpc3RhdGUgPT0gbnVsbCkgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIHJldHVybiBpc3RhdGUuaW5mbGF0ZVNldERpY3Rpb25hcnkodGhpcywgZGljdGlvbmFyeSwgZGljdExlbmd0aCk7XG59XG5cbi8qXG5cbiAgcHVibGljIGludCBkZWZsYXRlSW5pdChpbnQgbGV2ZWwpe1xuICAgIHJldHVybiBkZWZsYXRlSW5pdChsZXZlbCwgTUFYX1dCSVRTKTtcbiAgfVxuICBwdWJsaWMgaW50IGRlZmxhdGVJbml0KGludCBsZXZlbCwgYm9vbGVhbiBub3dyYXApe1xuICAgIHJldHVybiBkZWZsYXRlSW5pdChsZXZlbCwgTUFYX1dCSVRTLCBub3dyYXApO1xuICB9XG4gIHB1YmxpYyBpbnQgZGVmbGF0ZUluaXQoaW50IGxldmVsLCBpbnQgYml0cyl7XG4gICAgcmV0dXJuIGRlZmxhdGVJbml0KGxldmVsLCBiaXRzLCBmYWxzZSk7XG4gIH1cbiAgcHVibGljIGludCBkZWZsYXRlSW5pdChpbnQgbGV2ZWwsIGludCBiaXRzLCBib29sZWFuIG5vd3JhcCl7XG4gICAgZHN0YXRlPW5ldyBEZWZsYXRlKCk7XG4gICAgcmV0dXJuIGRzdGF0ZS5kZWZsYXRlSW5pdCh0aGlzLCBsZXZlbCwgbm93cmFwPy1iaXRzOmJpdHMpO1xuICB9XG4gIHB1YmxpYyBpbnQgZGVmbGF0ZShpbnQgZmx1c2gpe1xuICAgIGlmKGRzdGF0ZT09bnVsbCl7XG4gICAgICByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG4gICAgfVxuICAgIHJldHVybiBkc3RhdGUuZGVmbGF0ZSh0aGlzLCBmbHVzaCk7XG4gIH1cbiAgcHVibGljIGludCBkZWZsYXRlRW5kKCl7XG4gICAgaWYoZHN0YXRlPT1udWxsKSByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG4gICAgaW50IHJldD1kc3RhdGUuZGVmbGF0ZUVuZCgpO1xuICAgIGRzdGF0ZT1udWxsO1xuICAgIHJldHVybiByZXQ7XG4gIH1cbiAgcHVibGljIGludCBkZWZsYXRlUGFyYW1zKGludCBsZXZlbCwgaW50IHN0cmF0ZWd5KXtcbiAgICBpZihkc3RhdGU9PW51bGwpIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICByZXR1cm4gZHN0YXRlLmRlZmxhdGVQYXJhbXModGhpcywgbGV2ZWwsIHN0cmF0ZWd5KTtcbiAgfVxuICBwdWJsaWMgaW50IGRlZmxhdGVTZXREaWN0aW9uYXJ5IChieXRlW10gZGljdGlvbmFyeSwgaW50IGRpY3RMZW5ndGgpe1xuICAgIGlmKGRzdGF0ZSA9PSBudWxsKVxuICAgICAgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIHJldHVybiBkc3RhdGUuZGVmbGF0ZVNldERpY3Rpb25hcnkodGhpcywgZGljdGlvbmFyeSwgZGljdExlbmd0aCk7XG4gIH1cblxuKi9cblxuLypcbiAgLy8gRmx1c2ggYXMgbXVjaCBwZW5kaW5nIG91dHB1dCBhcyBwb3NzaWJsZS4gQWxsIGRlZmxhdGUoKSBvdXRwdXQgZ29lc1xuICAvLyB0aHJvdWdoIHRoaXMgZnVuY3Rpb24gc28gc29tZSBhcHBsaWNhdGlvbnMgbWF5IHdpc2ggdG8gbW9kaWZ5IGl0XG4gIC8vIHRvIGF2b2lkIGFsbG9jYXRpbmcgYSBsYXJnZSBzdHJtLT5uZXh0X291dCBidWZmZXIgYW5kIGNvcHlpbmcgaW50byBpdC5cbiAgLy8gKFNlZSBhbHNvIHJlYWRfYnVmKCkpLlxuICB2b2lkIGZsdXNoX3BlbmRpbmcoKXtcbiAgICBpbnQgbGVuPWRzdGF0ZS5wZW5kaW5nO1xuXG4gICAgaWYobGVuPmF2YWlsX291dCkgbGVuPWF2YWlsX291dDtcbiAgICBpZihsZW49PTApIHJldHVybjtcblxuICAgIGlmKGRzdGF0ZS5wZW5kaW5nX2J1Zi5sZW5ndGg8PWRzdGF0ZS5wZW5kaW5nX291dCB8fFxuICAgICAgIG5leHRfb3V0Lmxlbmd0aDw9bmV4dF9vdXRfaW5kZXggfHxcbiAgICAgICBkc3RhdGUucGVuZGluZ19idWYubGVuZ3RoPChkc3RhdGUucGVuZGluZ19vdXQrbGVuKSB8fFxuICAgICAgIG5leHRfb3V0Lmxlbmd0aDwobmV4dF9vdXRfaW5kZXgrbGVuKSl7XG4gICAgICBTeXN0ZW0ub3V0LnByaW50bG4oZHN0YXRlLnBlbmRpbmdfYnVmLmxlbmd0aCtcIiwgXCIrZHN0YXRlLnBlbmRpbmdfb3V0K1xuXHRcdFx0IFwiLCBcIituZXh0X291dC5sZW5ndGgrXCIsIFwiK25leHRfb3V0X2luZGV4K1wiLCBcIitsZW4pO1xuICAgICAgU3lzdGVtLm91dC5wcmludGxuKFwiYXZhaWxfb3V0PVwiK2F2YWlsX291dCk7XG4gICAgfVxuXG4gICAgU3lzdGVtLmFycmF5Y29weShkc3RhdGUucGVuZGluZ19idWYsIGRzdGF0ZS5wZW5kaW5nX291dCxcblx0XHQgICAgIG5leHRfb3V0LCBuZXh0X291dF9pbmRleCwgbGVuKTtcblxuICAgIG5leHRfb3V0X2luZGV4Kz1sZW47XG4gICAgZHN0YXRlLnBlbmRpbmdfb3V0Kz1sZW47XG4gICAgdG90YWxfb3V0Kz1sZW47XG4gICAgYXZhaWxfb3V0LT1sZW47XG4gICAgZHN0YXRlLnBlbmRpbmctPWxlbjtcbiAgICBpZihkc3RhdGUucGVuZGluZz09MCl7XG4gICAgICBkc3RhdGUucGVuZGluZ19vdXQ9MDtcbiAgICB9XG4gIH1cblxuICAvLyBSZWFkIGEgbmV3IGJ1ZmZlciBmcm9tIHRoZSBjdXJyZW50IGlucHV0IHN0cmVhbSwgdXBkYXRlIHRoZSBhZGxlcjMyXG4gIC8vIGFuZCB0b3RhbCBudW1iZXIgb2YgYnl0ZXMgcmVhZC4gIEFsbCBkZWZsYXRlKCkgaW5wdXQgZ29lcyB0aHJvdWdoXG4gIC8vIHRoaXMgZnVuY3Rpb24gc28gc29tZSBhcHBsaWNhdGlvbnMgbWF5IHdpc2ggdG8gbW9kaWZ5IGl0IHRvIGF2b2lkXG4gIC8vIGFsbG9jYXRpbmcgYSBsYXJnZSBzdHJtLT5uZXh0X2luIGJ1ZmZlciBhbmQgY29weWluZyBmcm9tIGl0LlxuICAvLyAoU2VlIGFsc28gZmx1c2hfcGVuZGluZygpKS5cbiAgaW50IHJlYWRfYnVmKGJ5dGVbXSBidWYsIGludCBzdGFydCwgaW50IHNpemUpIHtcbiAgICBpbnQgbGVuPWF2YWlsX2luO1xuXG4gICAgaWYobGVuPnNpemUpIGxlbj1zaXplO1xuICAgIGlmKGxlbj09MCkgcmV0dXJuIDA7XG5cbiAgICBhdmFpbF9pbi09bGVuO1xuXG4gICAgaWYoZHN0YXRlLm5vaGVhZGVyPT0wKSB7XG4gICAgICBhZGxlcj1fYWRsZXIuYWRsZXIzMihhZGxlciwgbmV4dF9pbiwgbmV4dF9pbl9pbmRleCwgbGVuKTtcbiAgICB9XG4gICAgU3lzdGVtLmFycmF5Y29weShuZXh0X2luLCBuZXh0X2luX2luZGV4LCBidWYsIHN0YXJ0LCBsZW4pO1xuICAgIG5leHRfaW5faW5kZXggICs9IGxlbjtcbiAgICB0b3RhbF9pbiArPSBsZW47XG4gICAgcmV0dXJuIGxlbjtcbiAgfVxuXG4gIHB1YmxpYyB2b2lkIGZyZWUoKXtcbiAgICBuZXh0X2luPW51bGw7XG4gICAgbmV4dF9vdXQ9bnVsbDtcbiAgICBtc2c9bnVsbDtcbiAgICBfYWRsZXI9bnVsbDtcbiAgfVxufVxuKi9cblxuXG4vL1xuLy8gSW5mbGF0ZS5qYXZhXG4vL1xuXG5mdW5jdGlvbiBJbmZsYXRlKCkge1xuICAgIHRoaXMud2FzID0gWzBdO1xufVxuXG5JbmZsYXRlLnByb3RvdHlwZS5pbmZsYXRlUmVzZXQgPSBmdW5jdGlvbih6KSB7XG4gICAgaWYoeiA9PSBudWxsIHx8IHouaXN0YXRlID09IG51bGwpIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICBcbiAgICB6LnRvdGFsX2luID0gei50b3RhbF9vdXQgPSAwO1xuICAgIHoubXNnID0gbnVsbDtcbiAgICB6LmlzdGF0ZS5tb2RlID0gei5pc3RhdGUubm93cmFwIT0wID8gQkxPQ0tTIDogTUVUSE9EO1xuICAgIHouaXN0YXRlLmJsb2Nrcy5yZXNldCh6LCBudWxsKTtcbiAgICByZXR1cm4gWl9PSztcbn1cblxuSW5mbGF0ZS5wcm90b3R5cGUuaW5mbGF0ZUVuZCA9IGZ1bmN0aW9uKHope1xuICAgIGlmKHRoaXMuYmxvY2tzICE9IG51bGwpXG4gICAgICB0aGlzLmJsb2Nrcy5mcmVlKHopO1xuICAgIHRoaXMuYmxvY2tzPW51bGw7XG4gICAgcmV0dXJuIFpfT0s7XG59XG5cbkluZmxhdGUucHJvdG90eXBlLmluZmxhdGVJbml0ID0gZnVuY3Rpb24oeiwgdyl7XG4gICAgei5tc2cgPSBudWxsO1xuICAgIHRoaXMuYmxvY2tzID0gbnVsbDtcblxuICAgIC8vIGhhbmRsZSB1bmRvY3VtZW50ZWQgbm93cmFwIG9wdGlvbiAobm8gemxpYiBoZWFkZXIgb3IgY2hlY2spXG4gICAgbm93cmFwID0gMDtcbiAgICBpZih3IDwgMCl7XG4gICAgICB3ID0gLSB3O1xuICAgICAgbm93cmFwID0gMTtcbiAgICB9XG5cbiAgICAvLyBzZXQgd2luZG93IHNpemVcbiAgICBpZih3PDggfHx3PjE1KXtcbiAgICAgIHRoaXMuaW5mbGF0ZUVuZCh6KTtcbiAgICAgIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICB9XG4gICAgdGhpcy53Yml0cz13O1xuXG4gICAgei5pc3RhdGUuYmxvY2tzPW5ldyBJbmZCbG9ja3MoeiwgXG5cdFx0XHRcdCAgei5pc3RhdGUubm93cmFwIT0wID8gbnVsbCA6IHRoaXMsXG5cdFx0XHRcdCAgMTw8dyk7XG5cbiAgICAvLyByZXNldCBzdGF0ZVxuICAgIHRoaXMuaW5mbGF0ZVJlc2V0KHopO1xuICAgIHJldHVybiBaX09LO1xuICB9XG5cbkluZmxhdGUucHJvdG90eXBlLmluZmxhdGUgPSBmdW5jdGlvbih6LCBmKXtcbiAgICB2YXIgciwgYjtcblxuICAgIGlmKHogPT0gbnVsbCB8fCB6LmlzdGF0ZSA9PSBudWxsIHx8IHoubmV4dF9pbiA9PSBudWxsKVxuICAgICAgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIGYgPSBmID09IFpfRklOSVNIID8gWl9CVUZfRVJST1IgOiBaX09LO1xuICAgIHIgPSBaX0JVRl9FUlJPUjtcbiAgICB3aGlsZSAodHJ1ZSl7XG4gICAgICBzd2l0Y2ggKHouaXN0YXRlLm1vZGUpe1xuICAgICAgY2FzZSBNRVRIT0Q6XG5cbiAgICAgICAgaWYoei5hdmFpbF9pbj09MClyZXR1cm4gcjtyPWY7XG5cbiAgICAgICAgei5hdmFpbF9pbi0tOyB6LnRvdGFsX2luKys7XG4gICAgICAgIGlmKCgoei5pc3RhdGUubWV0aG9kID0gei5uZXh0X2luW3oubmV4dF9pbl9pbmRleCsrXSkmMHhmKSE9Wl9ERUZMQVRFRCl7XG4gICAgICAgICAgei5pc3RhdGUubW9kZSA9IEJBRDtcbiAgICAgICAgICB6Lm1zZz1cInVua25vd24gY29tcHJlc3Npb24gbWV0aG9kXCI7XG4gICAgICAgICAgei5pc3RhdGUubWFya2VyID0gNTsgICAgICAgLy8gY2FuJ3QgdHJ5IGluZmxhdGVTeW5jXG4gICAgICAgICAgYnJlYWs7XG4gICAgICAgIH1cbiAgICAgICAgaWYoKHouaXN0YXRlLm1ldGhvZD4+NCkrOD56LmlzdGF0ZS53Yml0cyl7XG4gICAgICAgICAgei5pc3RhdGUubW9kZSA9IEJBRDtcbiAgICAgICAgICB6Lm1zZz1cImludmFsaWQgd2luZG93IHNpemVcIjtcbiAgICAgICAgICB6LmlzdGF0ZS5tYXJrZXIgPSA1OyAgICAgICAvLyBjYW4ndCB0cnkgaW5mbGF0ZVN5bmNcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgfVxuICAgICAgICB6LmlzdGF0ZS5tb2RlPUZMQUc7XG4gICAgICBjYXNlIEZMQUc6XG5cbiAgICAgICAgaWYoei5hdmFpbF9pbj09MClyZXR1cm4gcjtyPWY7XG5cbiAgICAgICAgei5hdmFpbF9pbi0tOyB6LnRvdGFsX2luKys7XG4gICAgICAgIGIgPSAoei5uZXh0X2luW3oubmV4dF9pbl9pbmRleCsrXSkmMHhmZjtcblxuICAgICAgICBpZigoKCh6LmlzdGF0ZS5tZXRob2QgPDwgOCkrYikgJSAzMSkhPTApe1xuICAgICAgICAgIHouaXN0YXRlLm1vZGUgPSBCQUQ7XG4gICAgICAgICAgei5tc2cgPSBcImluY29ycmVjdCBoZWFkZXIgY2hlY2tcIjtcbiAgICAgICAgICB6LmlzdGF0ZS5tYXJrZXIgPSA1OyAgICAgICAvLyBjYW4ndCB0cnkgaW5mbGF0ZVN5bmNcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgfVxuXG4gICAgICAgIGlmKChiJlBSRVNFVF9ESUNUKT09MCl7XG4gICAgICAgICAgei5pc3RhdGUubW9kZSA9IEJMT0NLUztcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgfVxuICAgICAgICB6LmlzdGF0ZS5tb2RlID0gRElDVDQ7XG4gICAgICBjYXNlIERJQ1Q0OlxuXG4gICAgICAgIGlmKHouYXZhaWxfaW49PTApcmV0dXJuIHI7cj1mO1xuXG4gICAgICAgIHouYXZhaWxfaW4tLTsgei50b3RhbF9pbisrO1xuICAgICAgICB6LmlzdGF0ZS5uZWVkPSgoei5uZXh0X2luW3oubmV4dF9pbl9pbmRleCsrXSYweGZmKTw8MjQpJjB4ZmYwMDAwMDA7XG4gICAgICAgIHouaXN0YXRlLm1vZGU9RElDVDM7XG4gICAgICBjYXNlIERJQ1QzOlxuXG4gICAgICAgIGlmKHouYXZhaWxfaW49PTApcmV0dXJuIHI7cj1mO1xuXG4gICAgICAgIHouYXZhaWxfaW4tLTsgei50b3RhbF9pbisrO1xuICAgICAgICB6LmlzdGF0ZS5uZWVkKz0oKHoubmV4dF9pblt6Lm5leHRfaW5faW5kZXgrK10mMHhmZik8PDE2KSYweGZmMDAwMDtcbiAgICAgICAgei5pc3RhdGUubW9kZT1ESUNUMjtcbiAgICAgIGNhc2UgRElDVDI6XG5cbiAgICAgICAgaWYoei5hdmFpbF9pbj09MClyZXR1cm4gcjtyPWY7XG5cbiAgICAgICAgei5hdmFpbF9pbi0tOyB6LnRvdGFsX2luKys7XG4gICAgICAgIHouaXN0YXRlLm5lZWQrPSgoei5uZXh0X2luW3oubmV4dF9pbl9pbmRleCsrXSYweGZmKTw8OCkmMHhmZjAwO1xuICAgICAgICB6LmlzdGF0ZS5tb2RlPURJQ1QxO1xuICAgICAgY2FzZSBESUNUMTpcblxuICAgICAgICBpZih6LmF2YWlsX2luPT0wKXJldHVybiByO3I9ZjtcblxuICAgICAgICB6LmF2YWlsX2luLS07IHoudG90YWxfaW4rKztcbiAgICAgICAgei5pc3RhdGUubmVlZCArPSAoei5uZXh0X2luW3oubmV4dF9pbl9pbmRleCsrXSYweGZmKTtcbiAgICAgICAgei5hZGxlciA9IHouaXN0YXRlLm5lZWQ7XG4gICAgICAgIHouaXN0YXRlLm1vZGUgPSBESUNUMDtcbiAgICAgICAgcmV0dXJuIFpfTkVFRF9ESUNUO1xuICAgICAgY2FzZSBESUNUMDpcbiAgICAgICAgei5pc3RhdGUubW9kZSA9IEJBRDtcbiAgICAgICAgei5tc2cgPSBcIm5lZWQgZGljdGlvbmFyeVwiO1xuICAgICAgICB6LmlzdGF0ZS5tYXJrZXIgPSAwOyAgICAgICAvLyBjYW4gdHJ5IGluZmxhdGVTeW5jXG4gICAgICAgIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICAgIGNhc2UgQkxPQ0tTOlxuXG4gICAgICAgIHIgPSB6LmlzdGF0ZS5ibG9ja3MucHJvYyh6LCByKTtcbiAgICAgICAgaWYociA9PSBaX0RBVEFfRVJST1Ipe1xuICAgICAgICAgIHouaXN0YXRlLm1vZGUgPSBCQUQ7XG4gICAgICAgICAgei5pc3RhdGUubWFya2VyID0gMDsgICAgIC8vIGNhbiB0cnkgaW5mbGF0ZVN5bmNcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgfVxuICAgICAgICBpZihyID09IFpfT0spe1xuICAgICAgICAgIHIgPSBmO1xuICAgICAgICB9XG4gICAgICAgIGlmKHIgIT0gWl9TVFJFQU1fRU5EKXtcbiAgICAgICAgICByZXR1cm4gcjtcbiAgICAgICAgfVxuICAgICAgICByID0gZjtcbiAgICAgICAgei5pc3RhdGUuYmxvY2tzLnJlc2V0KHosIHouaXN0YXRlLndhcyk7XG4gICAgICAgIGlmKHouaXN0YXRlLm5vd3JhcCE9MCl7XG4gICAgICAgICAgei5pc3RhdGUubW9kZT1ET05FO1xuICAgICAgICAgIGJyZWFrO1xuICAgICAgICB9XG4gICAgICAgIHouaXN0YXRlLm1vZGU9Q0hFQ0s0O1xuICAgICAgY2FzZSBDSEVDSzQ6XG5cbiAgICAgICAgaWYoei5hdmFpbF9pbj09MClyZXR1cm4gcjtyPWY7XG5cbiAgICAgICAgei5hdmFpbF9pbi0tOyB6LnRvdGFsX2luKys7XG4gICAgICAgIHouaXN0YXRlLm5lZWQ9KCh6Lm5leHRfaW5bei5uZXh0X2luX2luZGV4KytdJjB4ZmYpPDwyNCkmMHhmZjAwMDAwMDtcbiAgICAgICAgei5pc3RhdGUubW9kZT1DSEVDSzM7XG4gICAgICBjYXNlIENIRUNLMzpcblxuICAgICAgICBpZih6LmF2YWlsX2luPT0wKXJldHVybiByO3I9ZjtcblxuICAgICAgICB6LmF2YWlsX2luLS07IHoudG90YWxfaW4rKztcbiAgICAgICAgei5pc3RhdGUubmVlZCs9KCh6Lm5leHRfaW5bei5uZXh0X2luX2luZGV4KytdJjB4ZmYpPDwxNikmMHhmZjAwMDA7XG4gICAgICAgIHouaXN0YXRlLm1vZGUgPSBDSEVDSzI7XG4gICAgICBjYXNlIENIRUNLMjpcblxuICAgICAgICBpZih6LmF2YWlsX2luPT0wKXJldHVybiByO3I9ZjtcblxuICAgICAgICB6LmF2YWlsX2luLS07IHoudG90YWxfaW4rKztcbiAgICAgICAgei5pc3RhdGUubmVlZCs9KCh6Lm5leHRfaW5bei5uZXh0X2luX2luZGV4KytdJjB4ZmYpPDw4KSYweGZmMDA7XG4gICAgICAgIHouaXN0YXRlLm1vZGUgPSBDSEVDSzE7XG4gICAgICBjYXNlIENIRUNLMTpcblxuICAgICAgICBpZih6LmF2YWlsX2luPT0wKXJldHVybiByO3I9ZjtcblxuICAgICAgICB6LmF2YWlsX2luLS07IHoudG90YWxfaW4rKztcbiAgICAgICAgei5pc3RhdGUubmVlZCs9KHoubmV4dF9pblt6Lm5leHRfaW5faW5kZXgrK10mMHhmZik7XG5cbiAgICAgICAgaWYoKCh6LmlzdGF0ZS53YXNbMF0pKSAhPSAoKHouaXN0YXRlLm5lZWQpKSl7XG4gICAgICAgICAgei5pc3RhdGUubW9kZSA9IEJBRDtcbiAgICAgICAgICB6Lm1zZyA9IFwiaW5jb3JyZWN0IGRhdGEgY2hlY2tcIjtcbiAgICAgICAgICB6LmlzdGF0ZS5tYXJrZXIgPSA1OyAgICAgICAvLyBjYW4ndCB0cnkgaW5mbGF0ZVN5bmNcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgfVxuXG4gICAgICAgIHouaXN0YXRlLm1vZGUgPSBET05FO1xuICAgICAgY2FzZSBET05FOlxuICAgICAgICByZXR1cm4gWl9TVFJFQU1fRU5EO1xuICAgICAgY2FzZSBCQUQ6XG4gICAgICAgIHJldHVybiBaX0RBVEFfRVJST1I7XG4gICAgICBkZWZhdWx0OlxuICAgICAgICByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG4gICAgICB9XG4gICAgfVxuICB9XG5cblxuSW5mbGF0ZS5wcm90b3R5cGUuaW5mbGF0ZVNldERpY3Rpb25hcnkgPSBmdW5jdGlvbih6LCAgZGljdGlvbmFyeSwgZGljdExlbmd0aCkge1xuICAgIHZhciBpbmRleD0wO1xuICAgIHZhciBsZW5ndGggPSBkaWN0TGVuZ3RoO1xuICAgIGlmKHo9PW51bGwgfHwgei5pc3RhdGUgPT0gbnVsbHx8IHouaXN0YXRlLm1vZGUgIT0gRElDVDApXG4gICAgICByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG5cbiAgICBpZih6Ll9hZGxlci5hZGxlcjMyKDEsIGRpY3Rpb25hcnksIDAsIGRpY3RMZW5ndGgpIT16LmFkbGVyKXtcbiAgICAgIHJldHVybiBaX0RBVEFfRVJST1I7XG4gICAgfVxuXG4gICAgei5hZGxlciA9IHouX2FkbGVyLmFkbGVyMzIoMCwgbnVsbCwgMCwgMCk7XG5cbiAgICBpZihsZW5ndGggPj0gKDE8PHouaXN0YXRlLndiaXRzKSl7XG4gICAgICBsZW5ndGggPSAoMTw8ei5pc3RhdGUud2JpdHMpLTE7XG4gICAgICBpbmRleD1kaWN0TGVuZ3RoIC0gbGVuZ3RoO1xuICAgIH1cbiAgICB6LmlzdGF0ZS5ibG9ja3Muc2V0X2RpY3Rpb25hcnkoZGljdGlvbmFyeSwgaW5kZXgsIGxlbmd0aCk7XG4gICAgei5pc3RhdGUubW9kZSA9IEJMT0NLUztcbiAgICByZXR1cm4gWl9PSztcbiAgfVxuXG4vLyAgc3RhdGljIHByaXZhdGUgYnl0ZVtdIG1hcmsgPSB7KGJ5dGUpMCwgKGJ5dGUpMCwgKGJ5dGUpMHhmZiwgKGJ5dGUpMHhmZn07XG52YXIgbWFyayA9IFswLCAwLCAyNTUsIDI1NV1cblxuSW5mbGF0ZS5wcm90b3R5cGUuaW5mbGF0ZVN5bmMgPSBmdW5jdGlvbih6KXtcbiAgICB2YXIgbjsgICAgICAgLy8gbnVtYmVyIG9mIGJ5dGVzIHRvIGxvb2sgYXRcbiAgICB2YXIgcDsgICAgICAgLy8gcG9pbnRlciB0byBieXRlc1xuICAgIHZhciBtOyAgICAgICAvLyBudW1iZXIgb2YgbWFya2VyIGJ5dGVzIGZvdW5kIGluIGEgcm93XG4gICAgdmFyIHIsIHc7ICAgLy8gdGVtcG9yYXJpZXMgdG8gc2F2ZSB0b3RhbF9pbiBhbmQgdG90YWxfb3V0XG5cbiAgICAvLyBzZXQgdXBcbiAgICBpZih6ID09IG51bGwgfHwgei5pc3RhdGUgPT0gbnVsbClcbiAgICAgIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICBpZih6LmlzdGF0ZS5tb2RlICE9IEJBRCl7XG4gICAgICB6LmlzdGF0ZS5tb2RlID0gQkFEO1xuICAgICAgei5pc3RhdGUubWFya2VyID0gMDtcbiAgICB9XG4gICAgaWYoKG49ei5hdmFpbF9pbik9PTApXG4gICAgICByZXR1cm4gWl9CVUZfRVJST1I7XG4gICAgcD16Lm5leHRfaW5faW5kZXg7XG4gICAgbT16LmlzdGF0ZS5tYXJrZXI7XG5cbiAgICAvLyBzZWFyY2hcbiAgICB3aGlsZSAobiE9MCAmJiBtIDwgNCl7XG4gICAgICBpZih6Lm5leHRfaW5bcF0gPT0gbWFya1ttXSl7XG4gICAgICAgIG0rKztcbiAgICAgIH1cbiAgICAgIGVsc2UgaWYoei5uZXh0X2luW3BdIT0wKXtcbiAgICAgICAgbSA9IDA7XG4gICAgICB9XG4gICAgICBlbHNle1xuICAgICAgICBtID0gNCAtIG07XG4gICAgICB9XG4gICAgICBwKys7IG4tLTtcbiAgICB9XG5cbiAgICAvLyByZXN0b3JlXG4gICAgei50b3RhbF9pbiArPSBwLXoubmV4dF9pbl9pbmRleDtcbiAgICB6Lm5leHRfaW5faW5kZXggPSBwO1xuICAgIHouYXZhaWxfaW4gPSBuO1xuICAgIHouaXN0YXRlLm1hcmtlciA9IG07XG5cbiAgICAvLyByZXR1cm4gbm8gam95IG9yIHNldCB1cCB0byByZXN0YXJ0IG9uIGEgbmV3IGJsb2NrXG4gICAgaWYobSAhPSA0KXtcbiAgICAgIHJldHVybiBaX0RBVEFfRVJST1I7XG4gICAgfVxuICAgIHI9ei50b3RhbF9pbjsgIHc9ei50b3RhbF9vdXQ7XG4gICAgdGhpcy5pbmZsYXRlUmVzZXQoeik7XG4gICAgei50b3RhbF9pbj1yOyAgei50b3RhbF9vdXQgPSB3O1xuICAgIHouaXN0YXRlLm1vZGUgPSBCTE9DS1M7XG4gICAgcmV0dXJuIFpfT0s7XG59XG5cbiAgLy8gUmV0dXJucyB0cnVlIGlmIGluZmxhdGUgaXMgY3VycmVudGx5IGF0IHRoZSBlbmQgb2YgYSBibG9jayBnZW5lcmF0ZWRcbiAgLy8gYnkgWl9TWU5DX0ZMVVNIIG9yIFpfRlVMTF9GTFVTSC4gVGhpcyBmdW5jdGlvbiBpcyB1c2VkIGJ5IG9uZSBQUFBcbiAgLy8gaW1wbGVtZW50YXRpb24gdG8gcHJvdmlkZSBhbiBhZGRpdGlvbmFsIHNhZmV0eSBjaGVjay4gUFBQIHVzZXMgWl9TWU5DX0ZMVVNIXG4gIC8vIGJ1dCByZW1vdmVzIHRoZSBsZW5ndGggYnl0ZXMgb2YgdGhlIHJlc3VsdGluZyBlbXB0eSBzdG9yZWQgYmxvY2suIFdoZW5cbiAgLy8gZGVjb21wcmVzc2luZywgUFBQIGNoZWNrcyB0aGF0IGF0IHRoZSBlbmQgb2YgaW5wdXQgcGFja2V0LCBpbmZsYXRlIGlzXG4gIC8vIHdhaXRpbmcgZm9yIHRoZXNlIGxlbmd0aCBieXRlcy5cbkluZmxhdGUucHJvdG90eXBlLmluZmxhdGVTeW5jUG9pbnQgPSBmdW5jdGlvbih6KXtcbiAgICBpZih6ID09IG51bGwgfHwgei5pc3RhdGUgPT0gbnVsbCB8fCB6LmlzdGF0ZS5ibG9ja3MgPT0gbnVsbClcbiAgICAgIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICByZXR1cm4gei5pc3RhdGUuYmxvY2tzLnN5bmNfcG9pbnQoKTtcbn1cblxuXG4vL1xuLy8gSW5mQmxvY2tzLmphdmFcbi8vXG5cbnZhciBJTkZCTE9DS1NfQk9SREVSID0gWzE2LCAxNywgMTgsIDAsIDgsIDcsIDksIDYsIDEwLCA1LCAxMSwgNCwgMTIsIDMsIDEzLCAyLCAxNCwgMSwgMTVdO1xuXG5mdW5jdGlvbiBJbmZCbG9ja3MoeiwgY2hlY2tmbiwgdykge1xuICAgIHRoaXMuaHVmdHM9bmV3IEludDMyQXJyYXkoTUFOWSozKTtcbiAgICB0aGlzLndpbmRvdz1uZXcgVWludDhBcnJheSh3KTtcbiAgICB0aGlzLmVuZD13O1xuICAgIHRoaXMuY2hlY2tmbiA9IGNoZWNrZm47XG4gICAgdGhpcy5tb2RlID0gSUJfVFlQRTtcbiAgICB0aGlzLnJlc2V0KHosIG51bGwpO1xuXG4gICAgdGhpcy5sZWZ0ID0gMDsgICAgICAgICAgICAvLyBpZiBTVE9SRUQsIGJ5dGVzIGxlZnQgdG8gY29weSBcblxuICAgIHRoaXMudGFibGUgPSAwOyAgICAgICAgICAgLy8gdGFibGUgbGVuZ3RocyAoMTQgYml0cykgXG4gICAgdGhpcy5pbmRleCA9IDA7ICAgICAgICAgICAvLyBpbmRleCBpbnRvIGJsZW5zIChvciBib3JkZXIpIFxuICAgIHRoaXMuYmxlbnMgPSBudWxsOyAgICAgICAgIC8vIGJpdCBsZW5ndGhzIG9mIGNvZGVzIFxuICAgIHRoaXMuYmI9bmV3IEludDMyQXJyYXkoMSk7IC8vIGJpdCBsZW5ndGggdHJlZSBkZXB0aCBcbiAgICB0aGlzLnRiPW5ldyBJbnQzMkFycmF5KDEpOyAvLyBiaXQgbGVuZ3RoIGRlY29kaW5nIHRyZWUgXG5cbiAgICB0aGlzLmNvZGVzID0gbmV3IEluZkNvZGVzKCk7XG5cbiAgICB0aGlzLmxhc3QgPSAwOyAgICAgICAgICAgIC8vIHRydWUgaWYgdGhpcyBibG9jayBpcyB0aGUgbGFzdCBibG9jayBcblxuICAvLyBtb2RlIGluZGVwZW5kZW50IGluZm9ybWF0aW9uIFxuICAgIHRoaXMuYml0ayA9IDA7ICAgICAgICAgICAgLy8gYml0cyBpbiBiaXQgYnVmZmVyIFxuICAgIHRoaXMuYml0YiA9IDA7ICAgICAgICAgICAgLy8gYml0IGJ1ZmZlciBcbiAgICB0aGlzLnJlYWQgPSAwOyAgICAgICAgICAgIC8vIHdpbmRvdyByZWFkIHBvaW50ZXIgXG4gICAgdGhpcy53cml0ZSA9IDA7ICAgICAgICAgICAvLyB3aW5kb3cgd3JpdGUgcG9pbnRlciBcbiAgICB0aGlzLmNoZWNrID0gMDsgICAgICAgICAgLy8gY2hlY2sgb24gb3V0cHV0IFxuXG4gICAgdGhpcy5pbmZ0cmVlPW5ldyBJbmZUcmVlKCk7XG59XG5cblxuXG5cbkluZkJsb2Nrcy5wcm90b3R5cGUucmVzZXQgPSBmdW5jdGlvbih6LCBjKXtcbiAgICBpZihjKSBjWzBdPXRoaXMuY2hlY2s7XG4gICAgaWYodGhpcy5tb2RlPT1JQl9DT0RFUyl7XG4gICAgICB0aGlzLmNvZGVzLmZyZWUoeik7XG4gICAgfVxuICAgIHRoaXMubW9kZT1JQl9UWVBFO1xuICAgIHRoaXMuYml0az0wO1xuICAgIHRoaXMuYml0Yj0wO1xuICAgIHRoaXMucmVhZD10aGlzLndyaXRlPTA7XG5cbiAgICBpZih0aGlzLmNoZWNrZm4pXG4gICAgICB6LmFkbGVyPXRoaXMuY2hlY2s9ei5fYWRsZXIuYWRsZXIzMigwLCBudWxsLCAwLCAwKTtcbiAgfVxuXG4gSW5mQmxvY2tzLnByb3RvdHlwZS5wcm9jID0gZnVuY3Rpb24oeiwgcil7XG4gICAgdmFyIHQ7ICAgICAgICAgICAgICAvLyB0ZW1wb3Jhcnkgc3RvcmFnZVxuICAgIHZhciBiOyAgICAgICAgICAgICAgLy8gYml0IGJ1ZmZlclxuICAgIHZhciBrOyAgICAgICAgICAgICAgLy8gYml0cyBpbiBiaXQgYnVmZmVyXG4gICAgdmFyIHA7ICAgICAgICAgICAgICAvLyBpbnB1dCBkYXRhIHBvaW50ZXJcbiAgICB2YXIgbjsgICAgICAgICAgICAgIC8vIGJ5dGVzIGF2YWlsYWJsZSB0aGVyZVxuICAgIHZhciBxOyAgICAgICAgICAgICAgLy8gb3V0cHV0IHdpbmRvdyB3cml0ZSBwb2ludGVyXG4gICAgdmFyIG07ICAgICAgICAgICAgICAvLyBieXRlcyB0byBlbmQgb2Ygd2luZG93IG9yIHJlYWQgcG9pbnRlclxuXG4gICAgLy8gY29weSBpbnB1dC9vdXRwdXQgaW5mb3JtYXRpb24gdG8gbG9jYWxzIChVUERBVEUgbWFjcm8gcmVzdG9yZXMpXG4gICAge3A9ei5uZXh0X2luX2luZGV4O249ei5hdmFpbF9pbjtiPXRoaXMuYml0YjtrPXRoaXMuYml0azt9XG4gICAge3E9dGhpcy53cml0ZTttPShxPHRoaXMucmVhZCA/IHRoaXMucmVhZC1xLTEgOiB0aGlzLmVuZC1xKTt9XG5cbiAgICAvLyBwcm9jZXNzIGlucHV0IGJhc2VkIG9uIGN1cnJlbnQgc3RhdGVcbiAgICB3aGlsZSh0cnVlKXtcbiAgICAgIHN3aXRjaCAodGhpcy5tb2RlKXtcbiAgICAgIGNhc2UgSUJfVFlQRTpcblxuXHR3aGlsZShrPCgzKSl7XG5cdCAgaWYobiE9MCl7XG5cdCAgICByPVpfT0s7XG5cdCAgfVxuXHQgIGVsc2V7XG5cdCAgICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgICAgei5hdmFpbF9pbj1uO1xuXHQgICAgei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICB0aGlzLndyaXRlPXE7XG5cdCAgICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgfTtcblx0ICBuLS07XG5cdCAgYnw9KHoubmV4dF9pbltwKytdJjB4ZmYpPDxrO1xuXHQgIGsrPTg7XG5cdH1cblx0dCA9IChiICYgNyk7XG5cdHRoaXMubGFzdCA9IHQgJiAxO1xuXG5cdHN3aXRjaCAodCA+Pj4gMSl7XG4gICAgICAgIGNhc2UgMDogICAgICAgICAgICAgICAgICAgICAgICAgLy8gc3RvcmVkIFxuICAgICAgICAgIHtiPj4+PSgzKTtrLT0oMyk7fVxuICAgICAgICAgIHQgPSBrICYgNzsgICAgICAgICAgICAgICAgICAgIC8vIGdvIHRvIGJ5dGUgYm91bmRhcnlcblxuICAgICAgICAgIHtiPj4+PSh0KTtrLT0odCk7fVxuICAgICAgICAgIHRoaXMubW9kZSA9IElCX0xFTlM7ICAgICAgICAgICAgICAgICAgLy8gZ2V0IGxlbmd0aCBvZiBzdG9yZWQgYmxvY2tcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgY2FzZSAxOiAgICAgICAgICAgICAgICAgICAgICAgICAvLyBmaXhlZFxuICAgICAgICAgIHtcbiAgICAgICAgICAgICAgdmFyIGJsPW5ldyBJbnQzMkFycmF5KDEpO1xuXHQgICAgICB2YXIgYmQ9bmV3IEludDMyQXJyYXkoMSk7XG4gICAgICAgICAgICAgIHZhciB0bD1bXTtcblx0ICAgICAgdmFyIHRkPVtdO1xuXG5cdCAgICAgIGluZmxhdGVfdHJlZXNfZml4ZWQoYmwsIGJkLCB0bCwgdGQsIHopO1xuICAgICAgICAgICAgICB0aGlzLmNvZGVzLmluaXQoYmxbMF0sIGJkWzBdLCB0bFswXSwgMCwgdGRbMF0sIDAsIHopO1xuICAgICAgICAgIH1cblxuICAgICAgICAgIHtiPj4+PSgzKTtrLT0oMyk7fVxuXG4gICAgICAgICAgdGhpcy5tb2RlID0gSUJfQ09ERVM7XG4gICAgICAgICAgYnJlYWs7XG4gICAgICAgIGNhc2UgMjogICAgICAgICAgICAgICAgICAgICAgICAgLy8gZHluYW1pY1xuXG4gICAgICAgICAge2I+Pj49KDMpO2stPSgzKTt9XG5cbiAgICAgICAgICB0aGlzLm1vZGUgPSBJQl9UQUJMRTtcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgY2FzZSAzOiAgICAgICAgICAgICAgICAgICAgICAgICAvLyBpbGxlZ2FsXG5cbiAgICAgICAgICB7Yj4+Pj0oMyk7ay09KDMpO31cbiAgICAgICAgICB0aGlzLm1vZGUgPSBCQUQ7XG4gICAgICAgICAgei5tc2cgPSBcImludmFsaWQgYmxvY2sgdHlwZVwiO1xuICAgICAgICAgIHIgPSBaX0RBVEFfRVJST1I7XG5cblx0ICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICB0aGlzLndyaXRlPXE7XG5cdCAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHR9XG5cdGJyZWFrO1xuICAgICAgY2FzZSBJQl9MRU5TOlxuXHR3aGlsZShrPCgzMikpe1xuXHQgIGlmKG4hPTApe1xuXHQgICAgcj1aX09LO1xuXHQgIH1cblx0ICBlbHNle1xuXHQgICAgdGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0ICAgIHouYXZhaWxfaW49bjtcblx0ICAgIHoudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgICAgdGhpcy53cml0ZT1xO1xuXHQgICAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgIH07XG5cdCAgbi0tO1xuXHQgIGJ8PSh6Lm5leHRfaW5bcCsrXSYweGZmKTw8aztcblx0ICBrKz04O1xuXHR9XG5cblx0aWYgKCgoKH5iKSA+Pj4gMTYpICYgMHhmZmZmKSAhPSAoYiAmIDB4ZmZmZikpe1xuXHQgIHRoaXMubW9kZSA9IEJBRDtcblx0ICB6Lm1zZyA9IFwiaW52YWxpZCBzdG9yZWQgYmxvY2sgbGVuZ3Roc1wiO1xuXHQgIHIgPSBaX0RBVEFfRVJST1I7XG5cblx0ICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICB0aGlzLndyaXRlPXE7XG5cdCAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHR9XG5cdHRoaXMubGVmdCA9IChiICYgMHhmZmZmKTtcblx0YiA9IGsgPSAwOyAgICAgICAgICAgICAgICAgICAgICAgLy8gZHVtcCBiaXRzXG5cdHRoaXMubW9kZSA9IHRoaXMubGVmdCE9MCA/IElCX1NUT1JFRCA6ICh0aGlzLmxhc3QhPTAgPyBJQl9EUlkgOiBJQl9UWVBFKTtcblx0YnJlYWs7XG4gICAgICBjYXNlIElCX1NUT1JFRDpcblx0aWYgKG4gPT0gMCl7XG5cdCAgdGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0ICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgd3JpdGU9cTtcblx0ICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdH1cblxuXHRpZihtPT0wKXtcblx0ICBpZihxPT1lbmQmJnJlYWQhPTApe1xuXHQgICAgcT0wOyBtPShxPHRoaXMucmVhZCA/IHRoaXMucmVhZC1xLTEgOiB0aGlzLmVuZC1xKTtcblx0ICB9XG5cdCAgaWYobT09MCl7XG5cdCAgICB0aGlzLndyaXRlPXE7IFxuXHQgICAgcj10aGlzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICAgIHE9dGhpcy53cml0ZTsgbSA9IChxIDwgdGhpcy5yZWFkID8gdGhpcy5yZWFkLXEtMSA6IHRoaXMuZW5kLXEpO1xuXHQgICAgaWYocT09dGhpcy5lbmQgJiYgdGhpcy5yZWFkICE9IDApe1xuXHQgICAgICBxPTA7IG0gPSAocSA8IHRoaXMucmVhZCA/IHRoaXMucmVhZC1xLTEgOiB0aGlzLmVuZC1xKTtcblx0ICAgIH1cblx0ICAgIGlmKG09PTApe1xuXHQgICAgICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgICAgICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICAgIHRoaXMud3JpdGU9cTtcblx0ICAgICAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgfVxuXHQgIH1cblx0fVxuXHRyPVpfT0s7XG5cblx0dCA9IHRoaXMubGVmdDtcblx0aWYodD5uKSB0ID0gbjtcblx0aWYodD5tKSB0ID0gbTtcblx0YXJyYXlDb3B5KHoubmV4dF9pbiwgcCwgd2luZG93LCBxLCB0KTtcblx0cCArPSB0OyAgbiAtPSB0O1xuXHRxICs9IHQ7ICBtIC09IHQ7XG5cdGlmICgodGhpcy5sZWZ0IC09IHQpICE9IDApXG5cdCAgYnJlYWs7XG5cdHRoaXMubW9kZSA9ICh0aGlzLmxhc3QgIT0gMCA/IElCX0RSWSA6IElCX1RZUEUpO1xuXHRicmVhaztcbiAgICAgIGNhc2UgSUJfVEFCTEU6XG5cblx0d2hpbGUoazwoMTQpKXtcblx0ICBpZihuIT0wKXtcblx0ICAgIHI9Wl9PSztcblx0ICB9XG5cdCAgZWxzZXtcblx0ICAgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgICB6LmF2YWlsX2luPW47XG5cdCAgICB6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgIHRoaXMud3JpdGU9cTtcblx0ICAgIHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICB9O1xuXHQgIG4tLTtcblx0ICBifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7XG5cdCAgays9ODtcblx0fVxuXG5cdHRoaXMudGFibGUgPSB0ID0gKGIgJiAweDNmZmYpO1xuXHRpZiAoKHQgJiAweDFmKSA+IDI5IHx8ICgodCA+PiA1KSAmIDB4MWYpID4gMjkpXG5cdCAge1xuXHQgICAgdGhpcy5tb2RlID0gSUJfQkFEO1xuXHQgICAgei5tc2cgPSBcInRvbyBtYW55IGxlbmd0aCBvciBkaXN0YW5jZSBzeW1ib2xzXCI7XG5cdCAgICByID0gWl9EQVRBX0VSUk9SO1xuXG5cdCAgICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgICAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgICAgdGhpcy53cml0ZT1xO1xuXHQgICAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgIH1cblx0dCA9IDI1OCArICh0ICYgMHgxZikgKyAoKHQgPj4gNSkgJiAweDFmKTtcblx0aWYodGhpcy5ibGVucz09bnVsbCB8fCB0aGlzLmJsZW5zLmxlbmd0aDx0KXtcblx0ICAgIHRoaXMuYmxlbnM9bmV3IEludDMyQXJyYXkodCk7XG5cdH1cblx0ZWxzZXtcblx0ICBmb3IodmFyIGk9MDsgaTx0OyBpKyspe1xuICAgICAgICAgICAgICB0aGlzLmJsZW5zW2ldPTA7XG4gICAgICAgICAgfVxuXHR9XG5cblx0e2I+Pj49KDE0KTtrLT0oMTQpO31cblxuXHR0aGlzLmluZGV4ID0gMDtcblx0bW9kZSA9IElCX0JUUkVFO1xuICAgICAgY2FzZSBJQl9CVFJFRTpcblx0d2hpbGUgKHRoaXMuaW5kZXggPCA0ICsgKHRoaXMudGFibGUgPj4+IDEwKSl7XG5cdCAgd2hpbGUoazwoMykpe1xuXHQgICAgaWYobiE9MCl7XG5cdCAgICAgIHI9Wl9PSztcblx0ICAgIH1cblx0ICAgIGVsc2V7XG5cdCAgICAgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgICAgIHouYXZhaWxfaW49bjtcblx0ICAgICAgei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICAgIHRoaXMud3JpdGU9cTtcblx0ICAgICAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgfTtcblx0ICAgIG4tLTtcblx0ICAgIGJ8PSh6Lm5leHRfaW5bcCsrXSYweGZmKTw8aztcblx0ICAgIGsrPTg7XG5cdCAgfVxuXG5cdCAgdGhpcy5ibGVuc1tJTkZCTE9DS1NfQk9SREVSW3RoaXMuaW5kZXgrK11dID0gYiY3O1xuXG5cdCAge2I+Pj49KDMpO2stPSgzKTt9XG5cdH1cblxuXHR3aGlsZSh0aGlzLmluZGV4IDwgMTkpe1xuXHQgIHRoaXMuYmxlbnNbSU5GQkxPQ0tTX0JPUkRFUlt0aGlzLmluZGV4KytdXSA9IDA7XG5cdH1cblxuXHR0aGlzLmJiWzBdID0gNztcblx0dCA9IHRoaXMuaW5mdHJlZS5pbmZsYXRlX3RyZWVzX2JpdHModGhpcy5ibGVucywgdGhpcy5iYiwgdGhpcy50YiwgdGhpcy5odWZ0cywgeik7XG5cdGlmICh0ICE9IFpfT0spe1xuXHQgIHIgPSB0O1xuXHQgIGlmIChyID09IFpfREFUQV9FUlJPUil7XG5cdCAgICB0aGlzLmJsZW5zPW51bGw7XG5cdCAgICB0aGlzLm1vZGUgPSBJQl9CQUQ7XG5cdCAgfVxuXG5cdCAgdGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0ICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgd3JpdGU9cTtcblx0ICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdH1cblxuXHR0aGlzLmluZGV4ID0gMDtcblx0dGhpcy5tb2RlID0gSUJfRFRSRUU7XG4gICAgICBjYXNlIElCX0RUUkVFOlxuXHR3aGlsZSAodHJ1ZSl7XG5cdCAgdCA9IHRoaXMudGFibGU7XG5cdCAgaWYoISh0aGlzLmluZGV4IDwgMjU4ICsgKHQgJiAweDFmKSArICgodCA+PiA1KSAmIDB4MWYpKSl7XG5cdCAgICBicmVhaztcblx0ICB9XG5cblx0ICB2YXIgaDsgLy9pbnRbXVxuXHQgIHZhciBpLCBqLCBjO1xuXG5cdCAgdCA9IHRoaXMuYmJbMF07XG5cblx0ICB3aGlsZShrPCh0KSl7XG5cdCAgICBpZihuIT0wKXtcblx0ICAgICAgcj1aX09LO1xuXHQgICAgfVxuXHQgICAgZWxzZXtcblx0ICAgICAgdGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0ICAgICAgei5hdmFpbF9pbj1uO1xuXHQgICAgICB6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgICAgdGhpcy53cml0ZT1xO1xuXHQgICAgICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgICB9O1xuXHQgICAgbi0tO1xuXHQgICAgYnw9KHoubmV4dF9pbltwKytdJjB4ZmYpPDxrO1xuXHQgICAgays9ODtcblx0ICB9XG5cbi8vXHQgIGlmICh0aGlzLnRiWzBdPT0tMSl7XG4vLyAgICAgICAgICAgIGRsb2coXCJudWxsLi4uXCIpO1xuLy9cdCAgfVxuXG5cdCAgdD10aGlzLmh1ZnRzWyh0aGlzLnRiWzBdKyhiICYgaW5mbGF0ZV9tYXNrW3RdKSkqMysxXTtcblx0ICBjPXRoaXMuaHVmdHNbKHRoaXMudGJbMF0rKGIgJiBpbmZsYXRlX21hc2tbdF0pKSozKzJdO1xuXG5cdCAgaWYgKGMgPCAxNil7XG5cdCAgICBiPj4+PSh0KTtrLT0odCk7XG5cdCAgICB0aGlzLmJsZW5zW3RoaXMuaW5kZXgrK10gPSBjO1xuXHQgIH1cblx0ICBlbHNlIHsgLy8gYyA9PSAxNi4uMThcblx0ICAgIGkgPSBjID09IDE4ID8gNyA6IGMgLSAxNDtcblx0ICAgIGogPSBjID09IDE4ID8gMTEgOiAzO1xuXG5cdCAgICB3aGlsZShrPCh0K2kpKXtcblx0ICAgICAgaWYobiE9MCl7XG5cdFx0cj1aX09LO1xuXHQgICAgICB9XG5cdCAgICAgIGVsc2V7XG5cdFx0dGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0XHR6LmF2YWlsX2luPW47XG5cdFx0ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdFx0dGhpcy53cml0ZT1xO1xuXHRcdHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICAgICAgfTtcblx0ICAgICAgbi0tO1xuXHQgICAgICBifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7XG5cdCAgICAgIGsrPTg7XG5cdCAgICB9XG5cblx0ICAgIGI+Pj49KHQpO2stPSh0KTtcblxuXHQgICAgaiArPSAoYiAmIGluZmxhdGVfbWFza1tpXSk7XG5cblx0ICAgIGI+Pj49KGkpO2stPShpKTtcblxuXHQgICAgaSA9IHRoaXMuaW5kZXg7XG5cdCAgICB0ID0gdGhpcy50YWJsZTtcblx0ICAgIGlmIChpICsgaiA+IDI1OCArICh0ICYgMHgxZikgKyAoKHQgPj4gNSkgJiAweDFmKSB8fFxuXHRcdChjID09IDE2ICYmIGkgPCAxKSl7XG5cdCAgICAgIHRoaXMuYmxlbnM9bnVsbDtcblx0ICAgICAgdGhpcy5tb2RlID0gSUJfQkFEO1xuXHQgICAgICB6Lm1zZyA9IFwiaW52YWxpZCBiaXQgbGVuZ3RoIHJlcGVhdFwiO1xuXHQgICAgICByID0gWl9EQVRBX0VSUk9SO1xuXG5cdCAgICAgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgICAgdGhpcy53cml0ZT1xO1xuXHQgICAgICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgICB9XG5cblx0ICAgIGMgPSBjID09IDE2ID8gdGhpcy5ibGVuc1tpLTFdIDogMDtcblx0ICAgIGRve1xuXHQgICAgICB0aGlzLmJsZW5zW2krK10gPSBjO1xuXHQgICAgfVxuXHQgICAgd2hpbGUgKC0taiE9MCk7XG5cdCAgICB0aGlzLmluZGV4ID0gaTtcblx0ICB9XG5cdH1cblxuXHR0aGlzLnRiWzBdPS0xO1xuXHR7XG5cdCAgICB2YXIgYmw9bmV3IEludDMyQXJyYXkoMSk7XG5cdCAgICB2YXIgYmQ9bmV3IEludDMyQXJyYXkoMSk7XG5cdCAgICB2YXIgdGw9bmV3IEludDMyQXJyYXkoMSk7XG5cdCAgICB2YXIgdGQ9bmV3IEludDMyQXJyYXkoMSk7XG5cdCAgICBibFswXSA9IDk7ICAgICAgICAgLy8gbXVzdCBiZSA8PSA5IGZvciBsb29rYWhlYWQgYXNzdW1wdGlvbnNcblx0ICAgIGJkWzBdID0gNjsgICAgICAgICAvLyBtdXN0IGJlIDw9IDkgZm9yIGxvb2thaGVhZCBhc3N1bXB0aW9uc1xuXG5cdCAgICB0ID0gdGhpcy50YWJsZTtcblx0ICAgIHQgPSB0aGlzLmluZnRyZWUuaW5mbGF0ZV90cmVlc19keW5hbWljKDI1NyArICh0ICYgMHgxZiksIFxuXHRcdFx0XHRcdCAgICAgIDEgKyAoKHQgPj4gNSkgJiAweDFmKSxcblx0XHRcdFx0XHQgICAgICB0aGlzLmJsZW5zLCBibCwgYmQsIHRsLCB0ZCwgdGhpcy5odWZ0cywgeik7XG5cblx0ICAgIGlmICh0ICE9IFpfT0spe1xuXHQgICAgICAgIGlmICh0ID09IFpfREFUQV9FUlJPUil7XG5cdCAgICAgICAgICAgIHRoaXMuYmxlbnM9bnVsbDtcblx0ICAgICAgICAgICAgdGhpcy5tb2RlID0gQkFEO1xuXHQgICAgICAgIH1cblx0ICAgICAgICByID0gdDtcblxuXHQgICAgICAgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgICAgICAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgICAgICAgIHRoaXMud3JpdGU9cTtcblx0ICAgICAgICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgICB9XG5cdCAgICB0aGlzLmNvZGVzLmluaXQoYmxbMF0sIGJkWzBdLCB0aGlzLmh1ZnRzLCB0bFswXSwgdGhpcy5odWZ0cywgdGRbMF0sIHopO1xuXHR9XG5cdHRoaXMubW9kZSA9IElCX0NPREVTO1xuICAgICAgY2FzZSBJQl9DT0RFUzpcblx0dGhpcy5iaXRiPWI7IHRoaXMuYml0az1rO1xuXHR6LmF2YWlsX2luPW47IHoudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHR0aGlzLndyaXRlPXE7XG5cblx0aWYgKChyID0gdGhpcy5jb2Rlcy5wcm9jKHRoaXMsIHosIHIpKSAhPSBaX1NUUkVBTV9FTkQpe1xuXHQgIHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeiwgcik7XG5cdH1cblx0ciA9IFpfT0s7XG5cdHRoaXMuY29kZXMuZnJlZSh6KTtcblxuXHRwPXoubmV4dF9pbl9pbmRleDsgbj16LmF2YWlsX2luO2I9dGhpcy5iaXRiO2s9dGhpcy5iaXRrO1xuXHRxPXRoaXMud3JpdGU7bSA9IChxIDwgdGhpcy5yZWFkID8gdGhpcy5yZWFkLXEtMSA6IHRoaXMuZW5kLXEpO1xuXG5cdGlmICh0aGlzLmxhc3Q9PTApe1xuXHQgIHRoaXMubW9kZSA9IElCX1RZUEU7XG5cdCAgYnJlYWs7XG5cdH1cblx0dGhpcy5tb2RlID0gSUJfRFJZO1xuICAgICAgY2FzZSBJQl9EUlk6XG5cdHRoaXMud3JpdGU9cTsgXG5cdHIgPSB0aGlzLmluZmxhdGVfZmx1c2goeiwgcik7IFxuXHRxPXRoaXMud3JpdGU7IG0gPSAocSA8IHRoaXMucmVhZCA/IHRoaXMucmVhZC1xLTEgOiB0aGlzLmVuZC1xKTtcblx0aWYgKHRoaXMucmVhZCAhPSB0aGlzLndyaXRlKXtcblx0ICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICB0aGlzLndyaXRlPXE7XG5cdCAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LCByKTtcblx0fVxuXHRtb2RlID0gRE9ORTtcbiAgICAgIGNhc2UgSUJfRE9ORTpcblx0ciA9IFpfU1RSRUFNX0VORDtcblxuXHR0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHR6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdHRoaXMud3JpdGU9cTtcblx0cmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LCByKTtcbiAgICAgIGNhc2UgSUJfQkFEOlxuXHRyID0gWl9EQVRBX0VSUk9SO1xuXG5cdHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0dGhpcy53cml0ZT1xO1xuXHRyZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHosIHIpO1xuXG4gICAgICBkZWZhdWx0OlxuXHRyID0gWl9TVFJFQU1fRVJST1I7XG5cblx0dGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0ei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHR0aGlzLndyaXRlPXE7XG5cdHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeiwgcik7XG4gICAgICB9XG4gICAgfVxuICB9XG5cbkluZkJsb2Nrcy5wcm90b3R5cGUuZnJlZSA9IGZ1bmN0aW9uKHope1xuICAgIHRoaXMucmVzZXQoeiwgbnVsbCk7XG4gICAgdGhpcy53aW5kb3c9bnVsbDtcbiAgICB0aGlzLmh1ZnRzPW51bGw7XG59XG5cbkluZkJsb2Nrcy5wcm90b3R5cGUuc2V0X2RpY3Rpb25hcnkgPSBmdW5jdGlvbihkLCBzdGFydCwgbil7XG4gICAgYXJyYXlDb3B5KGQsIHN0YXJ0LCB3aW5kb3csIDAsIG4pO1xuICAgIHRoaXMucmVhZCA9IHRoaXMud3JpdGUgPSBuO1xufVxuXG4gIC8vIFJldHVybnMgdHJ1ZSBpZiBpbmZsYXRlIGlzIGN1cnJlbnRseSBhdCB0aGUgZW5kIG9mIGEgYmxvY2sgZ2VuZXJhdGVkXG4gIC8vIGJ5IFpfU1lOQ19GTFVTSCBvciBaX0ZVTExfRkxVU0guIFxuSW5mQmxvY2tzLnByb3RvdHlwZS5zeW5jX3BvaW50ID0gZnVuY3Rpb24oKXtcbiAgICByZXR1cm4gdGhpcy5tb2RlID09IElCX0xFTlM7XG59XG5cbiAgLy8gY29weSBhcyBtdWNoIGFzIHBvc3NpYmxlIGZyb20gdGhlIHNsaWRpbmcgd2luZG93IHRvIHRoZSBvdXRwdXQgYXJlYVxuSW5mQmxvY2tzLnByb3RvdHlwZS5pbmZsYXRlX2ZsdXNoID0gZnVuY3Rpb24oeiwgcil7XG4gICAgdmFyIG47XG4gICAgdmFyIHA7XG4gICAgdmFyIHE7XG5cbiAgICAvLyBsb2NhbCBjb3BpZXMgb2Ygc291cmNlIGFuZCBkZXN0aW5hdGlvbiBwb2ludGVyc1xuICAgIHAgPSB6Lm5leHRfb3V0X2luZGV4O1xuICAgIHEgPSB0aGlzLnJlYWQ7XG5cbiAgICAvLyBjb21wdXRlIG51bWJlciBvZiBieXRlcyB0byBjb3B5IGFzIGZhciBhcyBlbmQgb2Ygd2luZG93XG4gICAgbiA9ICgocSA8PSB0aGlzLndyaXRlID8gdGhpcy53cml0ZSA6IHRoaXMuZW5kKSAtIHEpO1xuICAgIGlmIChuID4gei5hdmFpbF9vdXQpIG4gPSB6LmF2YWlsX291dDtcbiAgICBpZiAobiE9MCAmJiByID09IFpfQlVGX0VSUk9SKSByID0gWl9PSztcblxuICAgIC8vIHVwZGF0ZSBjb3VudGVyc1xuICAgIHouYXZhaWxfb3V0IC09IG47XG4gICAgei50b3RhbF9vdXQgKz0gbjtcblxuICAgIC8vIHVwZGF0ZSBjaGVjayBpbmZvcm1hdGlvblxuICAgIGlmKHRoaXMuY2hlY2tmbiAhPSBudWxsKVxuICAgICAgei5hZGxlcj10aGlzLmNoZWNrPXouX2FkbGVyLmFkbGVyMzIodGhpcy5jaGVjaywgdGhpcy53aW5kb3csIHEsIG4pO1xuXG4gICAgLy8gY29weSBhcyBmYXIgYXMgZW5kIG9mIHdpbmRvd1xuICAgIGFycmF5Q29weSh0aGlzLndpbmRvdywgcSwgei5uZXh0X291dCwgcCwgbik7XG4gICAgcCArPSBuO1xuICAgIHEgKz0gbjtcblxuICAgIC8vIHNlZSBpZiBtb3JlIHRvIGNvcHkgYXQgYmVnaW5uaW5nIG9mIHdpbmRvd1xuICAgIGlmIChxID09IHRoaXMuZW5kKXtcbiAgICAgIC8vIHdyYXAgcG9pbnRlcnNcbiAgICAgIHEgPSAwO1xuICAgICAgaWYgKHRoaXMud3JpdGUgPT0gdGhpcy5lbmQpXG4gICAgICAgIHRoaXMud3JpdGUgPSAwO1xuXG4gICAgICAvLyBjb21wdXRlIGJ5dGVzIHRvIGNvcHlcbiAgICAgIG4gPSB0aGlzLndyaXRlIC0gcTtcbiAgICAgIGlmIChuID4gei5hdmFpbF9vdXQpIG4gPSB6LmF2YWlsX291dDtcbiAgICAgIGlmIChuIT0wICYmIHIgPT0gWl9CVUZfRVJST1IpIHIgPSBaX09LO1xuXG4gICAgICAvLyB1cGRhdGUgY291bnRlcnNcbiAgICAgIHouYXZhaWxfb3V0IC09IG47XG4gICAgICB6LnRvdGFsX291dCArPSBuO1xuXG4gICAgICAvLyB1cGRhdGUgY2hlY2sgaW5mb3JtYXRpb25cbiAgICAgIGlmKHRoaXMuY2hlY2tmbiAhPSBudWxsKVxuXHR6LmFkbGVyPXRoaXMuY2hlY2s9ei5fYWRsZXIuYWRsZXIzMih0aGlzLmNoZWNrLCB0aGlzLndpbmRvdywgcSwgbik7XG5cbiAgICAgIC8vIGNvcHlcbiAgICAgIGFycmF5Q29weSh0aGlzLndpbmRvdywgcSwgei5uZXh0X291dCwgcCwgbik7XG4gICAgICBwICs9IG47XG4gICAgICBxICs9IG47XG4gICAgfVxuXG4gICAgLy8gdXBkYXRlIHBvaW50ZXJzXG4gICAgei5uZXh0X291dF9pbmRleCA9IHA7XG4gICAgdGhpcy5yZWFkID0gcTtcblxuICAgIC8vIGRvbmVcbiAgICByZXR1cm4gcjtcbiAgfVxuXG4vL1xuLy8gSW5mQ29kZXMuamF2YVxuLy9cblxudmFyIElDX1NUQVJUPTA7ICAvLyB4OiBzZXQgdXAgZm9yIExFTlxudmFyIElDX0xFTj0xOyAgICAvLyBpOiBnZXQgbGVuZ3RoL2xpdGVyYWwvZW9iIG5leHRcbnZhciBJQ19MRU5FWFQ9MjsgLy8gaTogZ2V0dGluZyBsZW5ndGggZXh0cmEgKGhhdmUgYmFzZSlcbnZhciBJQ19ESVNUPTM7ICAgLy8gaTogZ2V0IGRpc3RhbmNlIG5leHRcbnZhciBJQ19ESVNURVhUPTQ7Ly8gaTogZ2V0dGluZyBkaXN0YW5jZSBleHRyYVxudmFyIElDX0NPUFk9NTsgICAvLyBvOiBjb3B5aW5nIGJ5dGVzIGluIHdpbmRvdywgd2FpdGluZyBmb3Igc3BhY2VcbnZhciBJQ19MSVQ9NjsgICAgLy8gbzogZ290IGxpdGVyYWwsIHdhaXRpbmcgZm9yIG91dHB1dCBzcGFjZVxudmFyIElDX1dBU0g9NzsgICAvLyBvOiBnb3QgZW9iLCBwb3NzaWJseSBzdGlsbCBvdXRwdXQgd2FpdGluZ1xudmFyIElDX0VORD04OyAgICAvLyB4OiBnb3QgZW9iIGFuZCBhbGwgZGF0YSBmbHVzaGVkXG52YXIgSUNfQkFEQ09ERT05Oy8vIHg6IGdvdCBlcnJvclxuXG5mdW5jdGlvbiBJbmZDb2RlcygpIHtcbn1cblxuSW5mQ29kZXMucHJvdG90eXBlLmluaXQgPSBmdW5jdGlvbihibCwgYmQsIHRsLCB0bF9pbmRleCwgdGQsIHRkX2luZGV4LCB6KSB7XG4gICAgdGhpcy5tb2RlPUlDX1NUQVJUO1xuICAgIHRoaXMubGJpdHM9Ymw7XG4gICAgdGhpcy5kYml0cz1iZDtcbiAgICB0aGlzLmx0cmVlPXRsO1xuICAgIHRoaXMubHRyZWVfaW5kZXg9dGxfaW5kZXg7XG4gICAgdGhpcy5kdHJlZSA9IHRkO1xuICAgIHRoaXMuZHRyZWVfaW5kZXg9dGRfaW5kZXg7XG4gICAgdGhpcy50cmVlPW51bGw7XG59XG5cbkluZkNvZGVzLnByb3RvdHlwZS5wcm9jID0gZnVuY3Rpb24ocywgeiwgcil7IFxuICAgIHZhciBqOyAgICAgICAgICAgICAgLy8gdGVtcG9yYXJ5IHN0b3JhZ2VcbiAgICB2YXIgdDsgICAgICAgICAgICAgIC8vIHRlbXBvcmFyeSBwb2ludGVyIChpbnRbXSlcbiAgICB2YXIgdGluZGV4OyAgICAgICAgIC8vIHRlbXBvcmFyeSBwb2ludGVyXG4gICAgdmFyIGU7ICAgICAgICAgICAgICAvLyBleHRyYSBiaXRzIG9yIG9wZXJhdGlvblxuICAgIHZhciBiPTA7ICAgICAgICAgICAgLy8gYml0IGJ1ZmZlclxuICAgIHZhciBrPTA7ICAgICAgICAgICAgLy8gYml0cyBpbiBiaXQgYnVmZmVyXG4gICAgdmFyIHA9MDsgICAgICAgICAgICAvLyBpbnB1dCBkYXRhIHBvaW50ZXJcbiAgICB2YXIgbjsgICAgICAgICAgICAgIC8vIGJ5dGVzIGF2YWlsYWJsZSB0aGVyZVxuICAgIHZhciBxOyAgICAgICAgICAgICAgLy8gb3V0cHV0IHdpbmRvdyB3cml0ZSBwb2ludGVyXG4gICAgdmFyIG07ICAgICAgICAgICAgICAvLyBieXRlcyB0byBlbmQgb2Ygd2luZG93IG9yIHJlYWQgcG9pbnRlclxuICAgIHZhciBmOyAgICAgICAgICAgICAgLy8gcG9pbnRlciB0byBjb3B5IHN0cmluZ3MgZnJvbVxuXG4gICAgLy8gY29weSBpbnB1dC9vdXRwdXQgaW5mb3JtYXRpb24gdG8gbG9jYWxzIChVUERBVEUgbWFjcm8gcmVzdG9yZXMpXG4gICAgcD16Lm5leHRfaW5faW5kZXg7bj16LmF2YWlsX2luO2I9cy5iaXRiO2s9cy5iaXRrO1xuICAgIHE9cy53cml0ZTttPXE8cy5yZWFkP3MucmVhZC1xLTE6cy5lbmQtcTtcblxuICAgIC8vIHByb2Nlc3MgaW5wdXQgYW5kIG91dHB1dCBiYXNlZCBvbiBjdXJyZW50IHN0YXRlXG4gICAgd2hpbGUgKHRydWUpe1xuICAgICAgc3dpdGNoICh0aGlzLm1vZGUpe1xuXHQvLyB3YWl0aW5nIGZvciBcImk6XCI9aW5wdXQsIFwibzpcIj1vdXRwdXQsIFwieDpcIj1ub3RoaW5nXG4gICAgICBjYXNlIElDX1NUQVJUOiAgICAgICAgIC8vIHg6IHNldCB1cCBmb3IgTEVOXG5cdGlmIChtID49IDI1OCAmJiBuID49IDEwKXtcblxuXHQgIHMuYml0Yj1iO3MuYml0az1rO1xuXHQgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICBzLndyaXRlPXE7XG5cdCAgciA9IHRoaXMuaW5mbGF0ZV9mYXN0KHRoaXMubGJpdHMsIHRoaXMuZGJpdHMsIFxuXHRcdFx0ICAgdGhpcy5sdHJlZSwgdGhpcy5sdHJlZV9pbmRleCwgXG5cdFx0XHQgICB0aGlzLmR0cmVlLCB0aGlzLmR0cmVlX2luZGV4LFxuXHRcdFx0ICAgcywgeik7XG5cblx0ICBwPXoubmV4dF9pbl9pbmRleDtuPXouYXZhaWxfaW47Yj1zLmJpdGI7az1zLmJpdGs7XG5cdCAgcT1zLndyaXRlO209cTxzLnJlYWQ/cy5yZWFkLXEtMTpzLmVuZC1xO1xuXG5cdCAgaWYgKHIgIT0gWl9PSyl7XG5cdCAgICB0aGlzLm1vZGUgPSByID09IFpfU1RSRUFNX0VORCA/IElDX1dBU0ggOiBJQ19CQURDT0RFO1xuXHQgICAgYnJlYWs7XG5cdCAgfVxuXHR9XG5cdHRoaXMubmVlZCA9IHRoaXMubGJpdHM7XG5cdHRoaXMudHJlZSA9IHRoaXMubHRyZWU7XG5cdHRoaXMudHJlZV9pbmRleD10aGlzLmx0cmVlX2luZGV4O1xuXG5cdHRoaXMubW9kZSA9IElDX0xFTjtcbiAgICAgIGNhc2UgSUNfTEVOOiAgICAgICAgICAgLy8gaTogZ2V0IGxlbmd0aC9saXRlcmFsL2VvYiBuZXh0XG5cdGogPSB0aGlzLm5lZWQ7XG5cblx0d2hpbGUoazwoaikpe1xuXHQgIGlmKG4hPTApcj1aX09LO1xuXHQgIGVsc2V7XG5cblx0ICAgIHMuYml0Yj1iO3MuYml0az1rO1xuXHQgICAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgICAgcy53cml0ZT1xO1xuXHQgICAgcmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgIH1cblx0ICBuLS07XG5cdCAgYnw9KHoubmV4dF9pbltwKytdJjB4ZmYpPDxrO1xuXHQgIGsrPTg7XG5cdH1cblxuXHR0aW5kZXg9KHRoaXMudHJlZV9pbmRleCsoYiZpbmZsYXRlX21hc2tbal0pKSozO1xuXG5cdGI+Pj49KHRoaXMudHJlZVt0aW5kZXgrMV0pO1xuXHRrLT0odGhpcy50cmVlW3RpbmRleCsxXSk7XG5cblx0ZT10aGlzLnRyZWVbdGluZGV4XTtcblxuXHRpZihlID09IDApeyAgICAgICAgICAgICAgIC8vIGxpdGVyYWxcblx0ICB0aGlzLmxpdCA9IHRoaXMudHJlZVt0aW5kZXgrMl07XG5cdCAgdGhpcy5tb2RlID0gSUNfTElUO1xuXHQgIGJyZWFrO1xuXHR9XG5cdGlmKChlICYgMTYpIT0wICl7ICAgICAgICAgIC8vIGxlbmd0aFxuXHQgIHRoaXMuZ2V0ID0gZSAmIDE1O1xuXHQgIHRoaXMubGVuID0gdGhpcy50cmVlW3RpbmRleCsyXTtcblx0ICB0aGlzLm1vZGUgPSBJQ19MRU5FWFQ7XG5cdCAgYnJlYWs7XG5cdH1cblx0aWYgKChlICYgNjQpID09IDApeyAgICAgICAgLy8gbmV4dCB0YWJsZVxuXHQgIHRoaXMubmVlZCA9IGU7XG5cdCAgdGhpcy50cmVlX2luZGV4ID0gdGluZGV4LzMgKyB0aGlzLnRyZWVbdGluZGV4KzJdO1xuXHQgIGJyZWFrO1xuXHR9XG5cdGlmICgoZSAmIDMyKSE9MCl7ICAgICAgICAgICAgICAgLy8gZW5kIG9mIGJsb2NrXG5cdCAgdGhpcy5tb2RlID0gSUNfV0FTSDtcblx0ICBicmVhaztcblx0fVxuXHR0aGlzLm1vZGUgPSBJQ19CQURDT0RFOyAgICAgICAgLy8gaW52YWxpZCBjb2RlXG5cdHoubXNnID0gXCJpbnZhbGlkIGxpdGVyYWwvbGVuZ3RoIGNvZGVcIjtcblx0ciA9IFpfREFUQV9FUlJPUjtcblxuXHRzLmJpdGI9YjtzLmJpdGs9aztcblx0ei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHRzLndyaXRlPXE7XG5cdHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcblxuICAgICAgY2FzZSBJQ19MRU5FWFQ6ICAgICAgICAvLyBpOiBnZXR0aW5nIGxlbmd0aCBleHRyYSAoaGF2ZSBiYXNlKVxuXHRqID0gdGhpcy5nZXQ7XG5cblx0d2hpbGUoazwoaikpe1xuXHQgIGlmKG4hPTApcj1aX09LO1xuXHQgIGVsc2V7XG5cblx0ICAgIHMuYml0Yj1iO3MuYml0az1rO1xuXHQgICAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgICAgcy53cml0ZT1xO1xuXHQgICAgcmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgIH1cblx0ICBuLS07IGJ8PSh6Lm5leHRfaW5bcCsrXSYweGZmKTw8aztcblx0ICBrKz04O1xuXHR9XG5cblx0dGhpcy5sZW4gKz0gKGIgJiBpbmZsYXRlX21hc2tbal0pO1xuXG5cdGI+Pj1qO1xuXHRrLT1qO1xuXG5cdHRoaXMubmVlZCA9IHRoaXMuZGJpdHM7XG5cdHRoaXMudHJlZSA9IHRoaXMuZHRyZWU7XG5cdHRoaXMudHJlZV9pbmRleCA9IHRoaXMuZHRyZWVfaW5kZXg7XG5cdHRoaXMubW9kZSA9IElDX0RJU1Q7XG4gICAgICBjYXNlIElDX0RJU1Q6ICAgICAgICAgIC8vIGk6IGdldCBkaXN0YW5jZSBuZXh0XG5cdGogPSB0aGlzLm5lZWQ7XG5cblx0d2hpbGUoazwoaikpe1xuXHQgIGlmKG4hPTApcj1aX09LO1xuXHQgIGVsc2V7XG5cblx0ICAgIHMuYml0Yj1iO3MuYml0az1rO1xuXHQgICAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgICAgcy53cml0ZT1xO1xuXHQgICAgcmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgIH1cblx0ICBuLS07IGJ8PSh6Lm5leHRfaW5bcCsrXSYweGZmKTw8aztcblx0ICBrKz04O1xuXHR9XG5cblx0dGluZGV4PSh0aGlzLnRyZWVfaW5kZXgrKGIgJiBpbmZsYXRlX21hc2tbal0pKSozO1xuXG5cdGI+Pj10aGlzLnRyZWVbdGluZGV4KzFdO1xuXHRrLT10aGlzLnRyZWVbdGluZGV4KzFdO1xuXG5cdGUgPSAodGhpcy50cmVlW3RpbmRleF0pO1xuXHRpZigoZSAmIDE2KSE9MCl7ICAgICAgICAgICAgICAgLy8gZGlzdGFuY2Vcblx0ICB0aGlzLmdldCA9IGUgJiAxNTtcblx0ICB0aGlzLmRpc3QgPSB0aGlzLnRyZWVbdGluZGV4KzJdO1xuXHQgIHRoaXMubW9kZSA9IElDX0RJU1RFWFQ7XG5cdCAgYnJlYWs7XG5cdH1cblx0aWYgKChlICYgNjQpID09IDApeyAgICAgICAgLy8gbmV4dCB0YWJsZVxuXHQgIHRoaXMubmVlZCA9IGU7XG5cdCAgdGhpcy50cmVlX2luZGV4ID0gdGluZGV4LzMgKyB0aGlzLnRyZWVbdGluZGV4KzJdO1xuXHQgIGJyZWFrO1xuXHR9XG5cdHRoaXMubW9kZSA9IElDX0JBRENPREU7ICAgICAgICAvLyBpbnZhbGlkIGNvZGVcblx0ei5tc2cgPSBcImludmFsaWQgZGlzdGFuY2UgY29kZVwiO1xuXHRyID0gWl9EQVRBX0VSUk9SO1xuXG5cdHMuYml0Yj1iO3MuYml0az1rO1xuXHR6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdHMud3JpdGU9cTtcblx0cmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXG4gICAgICBjYXNlIElDX0RJU1RFWFQ6ICAgICAgIC8vIGk6IGdldHRpbmcgZGlzdGFuY2UgZXh0cmFcblx0aiA9IHRoaXMuZ2V0O1xuXG5cdHdoaWxlKGs8KGopKXtcblx0ICBpZihuIT0wKXI9Wl9PSztcblx0ICBlbHNle1xuXG5cdCAgICBzLmJpdGI9YjtzLmJpdGs9aztcblx0ICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgIHMud3JpdGU9cTtcblx0ICAgIHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICB9XG5cdCAgbi0tOyBifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7XG5cdCAgays9ODtcblx0fVxuXG5cdHRoaXMuZGlzdCArPSAoYiAmIGluZmxhdGVfbWFza1tqXSk7XG5cblx0Yj4+PWo7XG5cdGstPWo7XG5cblx0dGhpcy5tb2RlID0gSUNfQ09QWTtcbiAgICAgIGNhc2UgSUNfQ09QWTogICAgICAgICAgLy8gbzogY29weWluZyBieXRlcyBpbiB3aW5kb3csIHdhaXRpbmcgZm9yIHNwYWNlXG4gICAgICAgIGYgPSBxIC0gdGhpcy5kaXN0O1xuICAgICAgICB3aGlsZShmIDwgMCl7ICAgICAvLyBtb2R1bG8gd2luZG93IHNpemUtXCJ3aGlsZVwiIGluc3RlYWRcbiAgICAgICAgICBmICs9IHMuZW5kOyAgICAgLy8gb2YgXCJpZlwiIGhhbmRsZXMgaW52YWxpZCBkaXN0YW5jZXNcblx0fVxuXHR3aGlsZSAodGhpcy5sZW4hPTApe1xuXG5cdCAgaWYobT09MCl7XG5cdCAgICBpZihxPT1zLmVuZCYmcy5yZWFkIT0wKXtxPTA7bT1xPHMucmVhZD9zLnJlYWQtcS0xOnMuZW5kLXE7fVxuXHQgICAgaWYobT09MCl7XG5cdCAgICAgIHMud3JpdGU9cTsgcj1zLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICAgICAgcT1zLndyaXRlO209cTxzLnJlYWQ/cy5yZWFkLXEtMTpzLmVuZC1xO1xuXG5cdCAgICAgIGlmKHE9PXMuZW5kJiZzLnJlYWQhPTApe3E9MDttPXE8cy5yZWFkP3MucmVhZC1xLTE6cy5lbmQtcTt9XG5cblx0ICAgICAgaWYobT09MCl7XG5cdFx0cy5iaXRiPWI7cy5iaXRrPWs7XG5cdFx0ei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHRcdHMud3JpdGU9cTtcblx0XHRyZXR1cm4gcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgICAgIH0gIFxuXHQgICAgfVxuXHQgIH1cblxuXHQgIHMud2luZG93W3ErK109cy53aW5kb3dbZisrXTsgbS0tO1xuXG5cdCAgaWYgKGYgPT0gcy5lbmQpXG4gICAgICAgICAgICBmID0gMDtcblx0ICB0aGlzLmxlbi0tO1xuXHR9XG5cdHRoaXMubW9kZSA9IElDX1NUQVJUO1xuXHRicmVhaztcbiAgICAgIGNhc2UgSUNfTElUOiAgICAgICAgICAgLy8gbzogZ290IGxpdGVyYWwsIHdhaXRpbmcgZm9yIG91dHB1dCBzcGFjZVxuXHRpZihtPT0wKXtcblx0ICBpZihxPT1zLmVuZCYmcy5yZWFkIT0wKXtxPTA7bT1xPHMucmVhZD9zLnJlYWQtcS0xOnMuZW5kLXE7fVxuXHQgIGlmKG09PTApe1xuXHQgICAgcy53cml0ZT1xOyByPXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgcT1zLndyaXRlO209cTxzLnJlYWQ/cy5yZWFkLXEtMTpzLmVuZC1xO1xuXG5cdCAgICBpZihxPT1zLmVuZCYmcy5yZWFkIT0wKXtxPTA7bT1xPHMucmVhZD9zLnJlYWQtcS0xOnMuZW5kLXE7fVxuXHQgICAgaWYobT09MCl7XG5cdCAgICAgIHMuYml0Yj1iO3MuYml0az1rO1xuXHQgICAgICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICAgIHMud3JpdGU9cTtcblx0ICAgICAgcmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgfVxuXHQgIH1cblx0fVxuXHRyPVpfT0s7XG5cblx0cy53aW5kb3dbcSsrXT10aGlzLmxpdDsgbS0tO1xuXG5cdHRoaXMubW9kZSA9IElDX1NUQVJUO1xuXHRicmVhaztcbiAgICAgIGNhc2UgSUNfV0FTSDogICAgICAgICAgIC8vIG86IGdvdCBlb2IsIHBvc3NpYmx5IG1vcmUgb3V0cHV0XG5cdGlmIChrID4gNyl7ICAgICAgICAvLyByZXR1cm4gdW51c2VkIGJ5dGUsIGlmIGFueVxuXHQgIGsgLT0gODtcblx0ICBuKys7XG5cdCAgcC0tOyAgICAgICAgICAgICAvLyBjYW4gYWx3YXlzIHJldHVybiBvbmVcblx0fVxuXG5cdHMud3JpdGU9cTsgcj1zLmluZmxhdGVfZmx1c2goeixyKTtcblx0cT1zLndyaXRlO209cTxzLnJlYWQ/cy5yZWFkLXEtMTpzLmVuZC1xO1xuXG5cdGlmIChzLnJlYWQgIT0gcy53cml0ZSl7XG5cdCAgcy5iaXRiPWI7cy5iaXRrPWs7XG5cdCAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgIHMud3JpdGU9cTtcblx0ICByZXR1cm4gcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdH1cblx0dGhpcy5tb2RlID0gSUNfRU5EO1xuICAgICAgY2FzZSBJQ19FTkQ6XG5cdHIgPSBaX1NUUkVBTV9FTkQ7XG5cdHMuYml0Yj1iO3MuYml0az1rO1xuXHR6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdHMud3JpdGU9cTtcblx0cmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXG4gICAgICBjYXNlIElDX0JBRENPREU6ICAgICAgIC8vIHg6IGdvdCBlcnJvclxuXG5cdHIgPSBaX0RBVEFfRVJST1I7XG5cblx0cy5iaXRiPWI7cy5iaXRrPWs7XG5cdHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0cy53cml0ZT1xO1xuXHRyZXR1cm4gcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cbiAgICAgIGRlZmF1bHQ6XG5cdHIgPSBaX1NUUkVBTV9FUlJPUjtcblxuXHRzLmJpdGI9YjtzLmJpdGs9aztcblx0ei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHRzLndyaXRlPXE7XG5cdHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcbiAgICAgIH1cbiAgICB9XG4gIH1cblxuSW5mQ29kZXMucHJvdG90eXBlLmZyZWUgPSBmdW5jdGlvbih6KXtcbiAgICAvLyAgWkZSRUUoeiwgYyk7XG59XG5cbiAgLy8gQ2FsbGVkIHdpdGggbnVtYmVyIG9mIGJ5dGVzIGxlZnQgdG8gd3JpdGUgaW4gd2luZG93IGF0IGxlYXN0IDI1OFxuICAvLyAodGhlIG1heGltdW0gc3RyaW5nIGxlbmd0aCkgYW5kIG51bWJlciBvZiBpbnB1dCBieXRlcyBhdmFpbGFibGVcbiAgLy8gYXQgbGVhc3QgdGVuLiAgVGhlIHRlbiBieXRlcyBhcmUgc2l4IGJ5dGVzIGZvciB0aGUgbG9uZ2VzdCBsZW5ndGgvXG4gIC8vIGRpc3RhbmNlIHBhaXIgcGx1cyBmb3VyIGJ5dGVzIGZvciBvdmVybG9hZGluZyB0aGUgYml0IGJ1ZmZlci5cblxuSW5mQ29kZXMucHJvdG90eXBlLmluZmxhdGVfZmFzdCA9IGZ1bmN0aW9uKGJsLCBiZCwgdGwsIHRsX2luZGV4LCB0ZCwgdGRfaW5kZXgsIHMsIHopIHtcbiAgICB2YXIgdDsgICAgICAgICAgICAgICAgLy8gdGVtcG9yYXJ5IHBvaW50ZXJcbiAgICB2YXIgICB0cDsgICAgICAgICAgICAgLy8gdGVtcG9yYXJ5IHBvaW50ZXIgKGludFtdKVxuICAgIHZhciB0cF9pbmRleDsgICAgICAgICAvLyB0ZW1wb3JhcnkgcG9pbnRlclxuICAgIHZhciBlOyAgICAgICAgICAgICAgICAvLyBleHRyYSBiaXRzIG9yIG9wZXJhdGlvblxuICAgIHZhciBiOyAgICAgICAgICAgICAgICAvLyBiaXQgYnVmZmVyXG4gICAgdmFyIGs7ICAgICAgICAgICAgICAgIC8vIGJpdHMgaW4gYml0IGJ1ZmZlclxuICAgIHZhciBwOyAgICAgICAgICAgICAgICAvLyBpbnB1dCBkYXRhIHBvaW50ZXJcbiAgICB2YXIgbjsgICAgICAgICAgICAgICAgLy8gYnl0ZXMgYXZhaWxhYmxlIHRoZXJlXG4gICAgdmFyIHE7ICAgICAgICAgICAgICAgIC8vIG91dHB1dCB3aW5kb3cgd3JpdGUgcG9pbnRlclxuICAgIHZhciBtOyAgICAgICAgICAgICAgICAvLyBieXRlcyB0byBlbmQgb2Ygd2luZG93IG9yIHJlYWQgcG9pbnRlclxuICAgIHZhciBtbDsgICAgICAgICAgICAgICAvLyBtYXNrIGZvciBsaXRlcmFsL2xlbmd0aCB0cmVlXG4gICAgdmFyIG1kOyAgICAgICAgICAgICAgIC8vIG1hc2sgZm9yIGRpc3RhbmNlIHRyZWVcbiAgICB2YXIgYzsgICAgICAgICAgICAgICAgLy8gYnl0ZXMgdG8gY29weVxuICAgIHZhciBkOyAgICAgICAgICAgICAgICAvLyBkaXN0YW5jZSBiYWNrIHRvIGNvcHkgZnJvbVxuICAgIHZhciByOyAgICAgICAgICAgICAgICAvLyBjb3B5IHNvdXJjZSBwb2ludGVyXG5cbiAgICB2YXIgdHBfaW5kZXhfdF8zOyAgICAgLy8gKHRwX2luZGV4K3QpKjNcblxuICAgIC8vIGxvYWQgaW5wdXQsIG91dHB1dCwgYml0IHZhbHVlc1xuICAgIHA9ei5uZXh0X2luX2luZGV4O249ei5hdmFpbF9pbjtiPXMuYml0YjtrPXMuYml0aztcbiAgICBxPXMud3JpdGU7bT1xPHMucmVhZD9zLnJlYWQtcS0xOnMuZW5kLXE7XG5cbiAgICAvLyBpbml0aWFsaXplIG1hc2tzXG4gICAgbWwgPSBpbmZsYXRlX21hc2tbYmxdO1xuICAgIG1kID0gaW5mbGF0ZV9tYXNrW2JkXTtcblxuICAgIC8vIGRvIHVudGlsIG5vdCBlbm91Z2ggaW5wdXQgb3Igb3V0cHV0IHNwYWNlIGZvciBmYXN0IGxvb3BcbiAgICBkbyB7ICAgICAgICAgICAgICAgICAgICAgICAgICAvLyBhc3N1bWUgY2FsbGVkIHdpdGggbSA+PSAyNTggJiYgbiA+PSAxMFxuICAgICAgLy8gZ2V0IGxpdGVyYWwvbGVuZ3RoIGNvZGVcbiAgICAgIHdoaWxlKGs8KDIwKSl7ICAgICAgICAgICAgICAvLyBtYXggYml0cyBmb3IgbGl0ZXJhbC9sZW5ndGggY29kZVxuXHRuLS07XG5cdGJ8PSh6Lm5leHRfaW5bcCsrXSYweGZmKTw8aztrKz04O1xuICAgICAgfVxuXG4gICAgICB0PSBiJm1sO1xuICAgICAgdHA9dGw7IFxuICAgICAgdHBfaW5kZXg9dGxfaW5kZXg7XG4gICAgICB0cF9pbmRleF90XzM9KHRwX2luZGV4K3QpKjM7XG4gICAgICBpZiAoKGUgPSB0cFt0cF9pbmRleF90XzNdKSA9PSAwKXtcblx0Yj4+PSh0cFt0cF9pbmRleF90XzMrMV0pOyBrLT0odHBbdHBfaW5kZXhfdF8zKzFdKTtcblxuXHRzLndpbmRvd1txKytdID0gdHBbdHBfaW5kZXhfdF8zKzJdO1xuXHRtLS07XG5cdGNvbnRpbnVlO1xuICAgICAgfVxuICAgICAgZG8ge1xuXG5cdGI+Pj0odHBbdHBfaW5kZXhfdF8zKzFdKTsgay09KHRwW3RwX2luZGV4X3RfMysxXSk7XG5cblx0aWYoKGUmMTYpIT0wKXtcblx0ICBlICY9IDE1O1xuXHQgIGMgPSB0cFt0cF9pbmRleF90XzMrMl0gKyAoYiAmIGluZmxhdGVfbWFza1tlXSk7XG5cblx0ICBiPj49ZTsgay09ZTtcblxuXHQgIC8vIGRlY29kZSBkaXN0YW5jZSBiYXNlIG9mIGJsb2NrIHRvIGNvcHlcblx0ICB3aGlsZShrPCgxNSkpeyAgICAgICAgICAgLy8gbWF4IGJpdHMgZm9yIGRpc3RhbmNlIGNvZGVcblx0ICAgIG4tLTtcblx0ICAgIGJ8PSh6Lm5leHRfaW5bcCsrXSYweGZmKTw8aztrKz04O1xuXHQgIH1cblxuXHQgIHQ9IGImbWQ7XG5cdCAgdHA9dGQ7XG5cdCAgdHBfaW5kZXg9dGRfaW5kZXg7XG4gICAgICAgICAgdHBfaW5kZXhfdF8zPSh0cF9pbmRleCt0KSozO1xuXHQgIGUgPSB0cFt0cF9pbmRleF90XzNdO1xuXG5cdCAgZG8ge1xuXG5cdCAgICBiPj49KHRwW3RwX2luZGV4X3RfMysxXSk7IGstPSh0cFt0cF9pbmRleF90XzMrMV0pO1xuXG5cdCAgICBpZigoZSYxNikhPTApe1xuXHQgICAgICAvLyBnZXQgZXh0cmEgYml0cyB0byBhZGQgdG8gZGlzdGFuY2UgYmFzZVxuXHQgICAgICBlICY9IDE1O1xuXHQgICAgICB3aGlsZShrPChlKSl7ICAgICAgICAgLy8gZ2V0IGV4dHJhIGJpdHMgKHVwIHRvIDEzKVxuXHRcdG4tLTtcblx0XHRifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7ays9ODtcblx0ICAgICAgfVxuXG5cdCAgICAgIGQgPSB0cFt0cF9pbmRleF90XzMrMl0gKyAoYiZpbmZsYXRlX21hc2tbZV0pO1xuXG5cdCAgICAgIGI+Pj0oZSk7IGstPShlKTtcblxuXHQgICAgICAvLyBkbyB0aGUgY29weVxuXHQgICAgICBtIC09IGM7XG5cdCAgICAgIGlmIChxID49IGQpeyAgICAgICAgICAgICAgICAvLyBvZmZzZXQgYmVmb3JlIGRlc3Rcblx0XHQvLyAganVzdCBjb3B5XG5cdFx0cj1xLWQ7XG5cdFx0aWYocS1yPjAgJiYgMj4ocS1yKSl7ICAgICAgICAgICBcblx0XHQgIHMud2luZG93W3ErK109cy53aW5kb3dbcisrXTsgLy8gbWluaW11bSBjb3VudCBpcyB0aHJlZSxcblx0XHQgIHMud2luZG93W3ErK109cy53aW5kb3dbcisrXTsgLy8gc28gdW5yb2xsIGxvb3AgYSBsaXR0bGVcblx0XHQgIGMtPTI7XG5cdFx0fVxuXHRcdGVsc2V7XG5cdFx0ICBzLndpbmRvd1txKytdPXMud2luZG93W3IrK107IC8vIG1pbmltdW0gY291bnQgaXMgdGhyZWUsXG5cdFx0ICBzLndpbmRvd1txKytdPXMud2luZG93W3IrK107IC8vIHNvIHVucm9sbCBsb29wIGEgbGl0dGxlXG5cdFx0ICBjLT0yO1xuXHRcdH1cblx0ICAgICAgfVxuXHQgICAgICBlbHNleyAgICAgICAgICAgICAgICAgIC8vIGVsc2Ugb2Zmc2V0IGFmdGVyIGRlc3RpbmF0aW9uXG4gICAgICAgICAgICAgICAgcj1xLWQ7XG4gICAgICAgICAgICAgICAgZG97XG4gICAgICAgICAgICAgICAgICByKz1zLmVuZDsgICAgICAgICAgLy8gZm9yY2UgcG9pbnRlciBpbiB3aW5kb3dcbiAgICAgICAgICAgICAgICB9d2hpbGUocjwwKTsgICAgICAgICAvLyBjb3ZlcnMgaW52YWxpZCBkaXN0YW5jZXNcblx0XHRlPXMuZW5kLXI7XG5cdFx0aWYoYz5lKXsgICAgICAgICAgICAgLy8gaWYgc291cmNlIGNyb3NzZXMsXG5cdFx0ICBjLT1lOyAgICAgICAgICAgICAgLy8gd3JhcHBlZCBjb3B5XG5cdFx0ICBpZihxLXI+MCAmJiBlPihxLXIpKXsgICAgICAgICAgIFxuXHRcdCAgICBkb3tzLndpbmRvd1txKytdID0gcy53aW5kb3dbcisrXTt9XG5cdFx0ICAgIHdoaWxlKC0tZSE9MCk7XG5cdFx0ICB9XG5cdFx0ICBlbHNle1xuXHRcdCAgICBhcnJheUNvcHkocy53aW5kb3csIHIsIHMud2luZG93LCBxLCBlKTtcblx0XHQgICAgcSs9ZTsgcis9ZTsgZT0wO1xuXHRcdCAgfVxuXHRcdCAgciA9IDA7ICAgICAgICAgICAgICAgICAgLy8gY29weSByZXN0IGZyb20gc3RhcnQgb2Ygd2luZG93XG5cdFx0fVxuXG5cdCAgICAgIH1cblxuXHQgICAgICAvLyBjb3B5IGFsbCBvciB3aGF0J3MgbGVmdFxuICAgICAgICAgICAgICBkb3tzLndpbmRvd1txKytdID0gcy53aW5kb3dbcisrXTt9XG5cdFx0d2hpbGUoLS1jIT0wKTtcblx0ICAgICAgYnJlYWs7XG5cdCAgICB9XG5cdCAgICBlbHNlIGlmKChlJjY0KT09MCl7XG5cdCAgICAgIHQrPXRwW3RwX2luZGV4X3RfMysyXTtcblx0ICAgICAgdCs9KGImaW5mbGF0ZV9tYXNrW2VdKTtcblx0ICAgICAgdHBfaW5kZXhfdF8zPSh0cF9pbmRleCt0KSozO1xuXHQgICAgICBlPXRwW3RwX2luZGV4X3RfM107XG5cdCAgICB9XG5cdCAgICBlbHNle1xuXHQgICAgICB6Lm1zZyA9IFwiaW52YWxpZCBkaXN0YW5jZSBjb2RlXCI7XG5cblx0ICAgICAgYz16LmF2YWlsX2luLW47Yz0oaz4+Myk8Yz9rPj4zOmM7bis9YztwLT1jO2stPWM8PDM7XG5cblx0ICAgICAgcy5iaXRiPWI7cy5iaXRrPWs7XG5cdCAgICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgICAgcy53cml0ZT1xO1xuXG5cdCAgICAgIHJldHVybiBaX0RBVEFfRVJST1I7XG5cdCAgICB9XG5cdCAgfVxuXHQgIHdoaWxlKHRydWUpO1xuXHQgIGJyZWFrO1xuXHR9XG5cblx0aWYoKGUmNjQpPT0wKXtcblx0ICB0Kz10cFt0cF9pbmRleF90XzMrMl07XG5cdCAgdCs9KGImaW5mbGF0ZV9tYXNrW2VdKTtcblx0ICB0cF9pbmRleF90XzM9KHRwX2luZGV4K3QpKjM7XG5cdCAgaWYoKGU9dHBbdHBfaW5kZXhfdF8zXSk9PTApe1xuXG5cdCAgICBiPj49KHRwW3RwX2luZGV4X3RfMysxXSk7IGstPSh0cFt0cF9pbmRleF90XzMrMV0pO1xuXG5cdCAgICBzLndpbmRvd1txKytdPXRwW3RwX2luZGV4X3RfMysyXTtcblx0ICAgIG0tLTtcblx0ICAgIGJyZWFrO1xuXHQgIH1cblx0fVxuXHRlbHNlIGlmKChlJjMyKSE9MCl7XG5cblx0ICBjPXouYXZhaWxfaW4tbjtjPShrPj4zKTxjP2s+PjM6YztuKz1jO3AtPWM7ay09Yzw8MztcbiBcblx0ICBzLmJpdGI9YjtzLmJpdGs9aztcblx0ICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgcy53cml0ZT1xO1xuXG5cdCAgcmV0dXJuIFpfU1RSRUFNX0VORDtcblx0fVxuXHRlbHNle1xuXHQgIHoubXNnPVwiaW52YWxpZCBsaXRlcmFsL2xlbmd0aCBjb2RlXCI7XG5cblx0ICBjPXouYXZhaWxfaW4tbjtjPShrPj4zKTxjP2s+PjM6YztuKz1jO3AtPWM7ay09Yzw8MztcblxuXHQgIHMuYml0Yj1iO3MuYml0az1rO1xuXHQgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICBzLndyaXRlPXE7XG5cblx0ICByZXR1cm4gWl9EQVRBX0VSUk9SO1xuXHR9XG4gICAgICB9IFxuICAgICAgd2hpbGUodHJ1ZSk7XG4gICAgfSBcbiAgICB3aGlsZShtPj0yNTggJiYgbj49IDEwKTtcblxuICAgIC8vIG5vdCBlbm91Z2ggaW5wdXQgb3Igb3V0cHV0LS1yZXN0b3JlIHBvaW50ZXJzIGFuZCByZXR1cm5cbiAgICBjPXouYXZhaWxfaW4tbjtjPShrPj4zKTxjP2s+PjM6YztuKz1jO3AtPWM7ay09Yzw8MztcblxuICAgIHMuYml0Yj1iO3MuYml0az1rO1xuICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcbiAgICBzLndyaXRlPXE7XG5cbiAgICByZXR1cm4gWl9PSztcbn1cblxuLy9cbi8vIEluZlRyZWUuamF2YVxuLy9cblxuZnVuY3Rpb24gSW5mVHJlZSgpIHtcbn1cblxuSW5mVHJlZS5wcm90b3R5cGUuaHVmdF9idWlsZCA9IGZ1bmN0aW9uKGIsIGJpbmRleCwgbiwgcywgZCwgZSwgdCwgbSwgaHAsIGhuLCB2KSB7XG5cbiAgICAvLyBHaXZlbiBhIGxpc3Qgb2YgY29kZSBsZW5ndGhzIGFuZCBhIG1heGltdW0gdGFibGUgc2l6ZSwgbWFrZSBhIHNldCBvZlxuICAgIC8vIHRhYmxlcyB0byBkZWNvZGUgdGhhdCBzZXQgb2YgY29kZXMuICBSZXR1cm4gWl9PSyBvbiBzdWNjZXNzLCBaX0JVRl9FUlJPUlxuICAgIC8vIGlmIHRoZSBnaXZlbiBjb2RlIHNldCBpcyBpbmNvbXBsZXRlICh0aGUgdGFibGVzIGFyZSBzdGlsbCBidWlsdCBpbiB0aGlzXG4gICAgLy8gY2FzZSksIFpfREFUQV9FUlJPUiBpZiB0aGUgaW5wdXQgaXMgaW52YWxpZCAoYW4gb3Zlci1zdWJzY3JpYmVkIHNldCBvZlxuICAgIC8vIGxlbmd0aHMpLCBvciBaX01FTV9FUlJPUiBpZiBub3QgZW5vdWdoIG1lbW9yeS5cblxuICAgIHZhciBhOyAgICAgICAgICAgICAgICAgICAgICAgLy8gY291bnRlciBmb3IgY29kZXMgb2YgbGVuZ3RoIGtcbiAgICB2YXIgZjsgICAgICAgICAgICAgICAgICAgICAgIC8vIGkgcmVwZWF0cyBpbiB0YWJsZSBldmVyeSBmIGVudHJpZXNcbiAgICB2YXIgZzsgICAgICAgICAgICAgICAgICAgICAgIC8vIG1heGltdW0gY29kZSBsZW5ndGhcbiAgICB2YXIgaDsgICAgICAgICAgICAgICAgICAgICAgIC8vIHRhYmxlIGxldmVsXG4gICAgdmFyIGk7ICAgICAgICAgICAgICAgICAgICAgICAvLyBjb3VudGVyLCBjdXJyZW50IGNvZGVcbiAgICB2YXIgajsgICAgICAgICAgICAgICAgICAgICAgIC8vIGNvdW50ZXJcbiAgICB2YXIgazsgICAgICAgICAgICAgICAgICAgICAgIC8vIG51bWJlciBvZiBiaXRzIGluIGN1cnJlbnQgY29kZVxuICAgIHZhciBsOyAgICAgICAgICAgICAgICAgICAgICAgLy8gYml0cyBwZXIgdGFibGUgKHJldHVybmVkIGluIG0pXG4gICAgdmFyIG1hc2s7ICAgICAgICAgICAgICAgICAgICAvLyAoMSA8PCB3KSAtIDEsIHRvIGF2b2lkIGNjIC1PIGJ1ZyBvbiBIUFxuICAgIHZhciBwOyAgICAgICAgICAgICAgICAgICAgICAgLy8gcG9pbnRlciBpbnRvIGNbXSwgYltdLCBvciB2W11cbiAgICB2YXIgcTsgICAgICAgICAgICAgICAgICAgICAgIC8vIHBvaW50cyB0byBjdXJyZW50IHRhYmxlXG4gICAgdmFyIHc7ICAgICAgICAgICAgICAgICAgICAgICAvLyBiaXRzIGJlZm9yZSB0aGlzIHRhYmxlID09IChsICogaClcbiAgICB2YXIgeHA7ICAgICAgICAgICAgICAgICAgICAgIC8vIHBvaW50ZXIgaW50byB4XG4gICAgdmFyIHk7ICAgICAgICAgICAgICAgICAgICAgICAvLyBudW1iZXIgb2YgZHVtbXkgY29kZXMgYWRkZWRcbiAgICB2YXIgejsgICAgICAgICAgICAgICAgICAgICAgIC8vIG51bWJlciBvZiBlbnRyaWVzIGluIGN1cnJlbnQgdGFibGVcblxuICAgIC8vIEdlbmVyYXRlIGNvdW50cyBmb3IgZWFjaCBiaXQgbGVuZ3RoXG5cbiAgICBwID0gMDsgaSA9IG47XG4gICAgZG8ge1xuICAgICAgdGhpcy5jW2JbYmluZGV4K3BdXSsrOyBwKys7IGktLTsgICAvLyBhc3N1bWUgYWxsIGVudHJpZXMgPD0gQk1BWFxuICAgIH13aGlsZShpIT0wKTtcblxuICAgIGlmKHRoaXMuY1swXSA9PSBuKXsgICAgICAgICAgICAgICAgLy8gbnVsbCBpbnB1dC0tYWxsIHplcm8gbGVuZ3RoIGNvZGVzXG4gICAgICB0WzBdID0gLTE7XG4gICAgICBtWzBdID0gMDtcbiAgICAgIHJldHVybiBaX09LO1xuICAgIH1cblxuICAgIC8vIEZpbmQgbWluaW11bSBhbmQgbWF4aW11bSBsZW5ndGgsIGJvdW5kICptIGJ5IHRob3NlXG4gICAgbCA9IG1bMF07XG4gICAgZm9yIChqID0gMTsgaiA8PSBCTUFYOyBqKyspXG4gICAgICBpZih0aGlzLmNbal0hPTApIGJyZWFrO1xuICAgIGsgPSBqOyAgICAgICAgICAgICAgICAgICAgICAgIC8vIG1pbmltdW0gY29kZSBsZW5ndGhcbiAgICBpZihsIDwgail7XG4gICAgICBsID0gajtcbiAgICB9XG4gICAgZm9yIChpID0gQk1BWDsgaSE9MDsgaS0tKXtcbiAgICAgIGlmKHRoaXMuY1tpXSE9MCkgYnJlYWs7XG4gICAgfVxuICAgIGcgPSBpOyAgICAgICAgICAgICAgICAgICAgICAgIC8vIG1heGltdW0gY29kZSBsZW5ndGhcbiAgICBpZihsID4gaSl7XG4gICAgICBsID0gaTtcbiAgICB9XG4gICAgbVswXSA9IGw7XG5cbiAgICAvLyBBZGp1c3QgbGFzdCBsZW5ndGggY291bnQgdG8gZmlsbCBvdXQgY29kZXMsIGlmIG5lZWRlZFxuICAgIGZvciAoeSA9IDEgPDwgajsgaiA8IGk7IGorKywgeSA8PD0gMSl7XG4gICAgICBpZiAoKHkgLT0gdGhpcy5jW2pdKSA8IDApe1xuICAgICAgICByZXR1cm4gWl9EQVRBX0VSUk9SO1xuICAgICAgfVxuICAgIH1cbiAgICBpZiAoKHkgLT0gdGhpcy5jW2ldKSA8IDApe1xuICAgICAgcmV0dXJuIFpfREFUQV9FUlJPUjtcbiAgICB9XG4gICAgdGhpcy5jW2ldICs9IHk7XG5cbiAgICAvLyBHZW5lcmF0ZSBzdGFydGluZyBvZmZzZXRzIGludG8gdGhlIHZhbHVlIHRhYmxlIGZvciBlYWNoIGxlbmd0aFxuICAgIHRoaXMueFsxXSA9IGogPSAwO1xuICAgIHAgPSAxOyAgeHAgPSAyO1xuICAgIHdoaWxlICgtLWkhPTApIHsgICAgICAgICAgICAgICAgIC8vIG5vdGUgdGhhdCBpID09IGcgZnJvbSBhYm92ZVxuICAgICAgdGhpcy54W3hwXSA9IChqICs9IHRoaXMuY1twXSk7XG4gICAgICB4cCsrO1xuICAgICAgcCsrO1xuICAgIH1cblxuICAgIC8vIE1ha2UgYSB0YWJsZSBvZiB2YWx1ZXMgaW4gb3JkZXIgb2YgYml0IGxlbmd0aHNcbiAgICBpID0gMDsgcCA9IDA7XG4gICAgZG8ge1xuICAgICAgaWYgKChqID0gYltiaW5kZXgrcF0pICE9IDApe1xuICAgICAgICB0aGlzLnZbdGhpcy54W2pdKytdID0gaTtcbiAgICAgIH1cbiAgICAgIHArKztcbiAgICB9XG4gICAgd2hpbGUgKCsraSA8IG4pO1xuICAgIG4gPSB0aGlzLnhbZ107ICAgICAgICAgICAgICAgICAgICAgLy8gc2V0IG4gdG8gbGVuZ3RoIG9mIHZcblxuICAgIC8vIEdlbmVyYXRlIHRoZSBIdWZmbWFuIGNvZGVzIGFuZCBmb3IgZWFjaCwgbWFrZSB0aGUgdGFibGUgZW50cmllc1xuICAgIHRoaXMueFswXSA9IGkgPSAwOyAgICAgICAgICAgICAgICAgLy8gZmlyc3QgSHVmZm1hbiBjb2RlIGlzIHplcm9cbiAgICBwID0gMDsgICAgICAgICAgICAgICAgICAgICAgICAvLyBncmFiIHZhbHVlcyBpbiBiaXQgb3JkZXJcbiAgICBoID0gLTE7ICAgICAgICAgICAgICAgICAgICAgICAvLyBubyB0YWJsZXMgeWV0LS1sZXZlbCAtMVxuICAgIHcgPSAtbDsgICAgICAgICAgICAgICAgICAgICAgIC8vIGJpdHMgZGVjb2RlZCA9PSAobCAqIGgpXG4gICAgdGhpcy51WzBdID0gMDsgICAgICAgICAgICAgICAgICAgICAvLyBqdXN0IHRvIGtlZXAgY29tcGlsZXJzIGhhcHB5XG4gICAgcSA9IDA7ICAgICAgICAgICAgICAgICAgICAgICAgLy8gZGl0dG9cbiAgICB6ID0gMDsgICAgICAgICAgICAgICAgICAgICAgICAvLyBkaXR0b1xuXG4gICAgLy8gZ28gdGhyb3VnaCB0aGUgYml0IGxlbmd0aHMgKGsgYWxyZWFkeSBpcyBiaXRzIGluIHNob3J0ZXN0IGNvZGUpXG4gICAgZm9yICg7IGsgPD0gZzsgaysrKXtcbiAgICAgIGEgPSB0aGlzLmNba107XG4gICAgICB3aGlsZSAoYS0tIT0wKXtcblx0Ly8gaGVyZSBpIGlzIHRoZSBIdWZmbWFuIGNvZGUgb2YgbGVuZ3RoIGsgYml0cyBmb3IgdmFsdWUgKnBcblx0Ly8gbWFrZSB0YWJsZXMgdXAgdG8gcmVxdWlyZWQgbGV2ZWxcbiAgICAgICAgd2hpbGUgKGsgPiB3ICsgbCl7XG4gICAgICAgICAgaCsrO1xuICAgICAgICAgIHcgKz0gbDsgICAgICAgICAgICAgICAgIC8vIHByZXZpb3VzIHRhYmxlIGFsd2F5cyBsIGJpdHNcblx0ICAvLyBjb21wdXRlIG1pbmltdW0gc2l6ZSB0YWJsZSBsZXNzIHRoYW4gb3IgZXF1YWwgdG8gbCBiaXRzXG4gICAgICAgICAgeiA9IGcgLSB3O1xuICAgICAgICAgIHogPSAoeiA+IGwpID8gbCA6IHo7ICAgICAgICAvLyB0YWJsZSBzaXplIHVwcGVyIGxpbWl0XG4gICAgICAgICAgaWYoKGY9MTw8KGo9ay13KSk+YSsxKXsgICAgIC8vIHRyeSBhIGstdyBiaXQgdGFibGVcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgLy8gdG9vIGZldyBjb2RlcyBmb3Igay13IGJpdCB0YWJsZVxuICAgICAgICAgICAgZiAtPSBhICsgMTsgICAgICAgICAgICAgICAvLyBkZWR1Y3QgY29kZXMgZnJvbSBwYXR0ZXJucyBsZWZ0XG4gICAgICAgICAgICB4cCA9IGs7XG4gICAgICAgICAgICBpZihqIDwgeil7XG4gICAgICAgICAgICAgIHdoaWxlICgrK2ogPCB6KXsgICAgICAgIC8vIHRyeSBzbWFsbGVyIHRhYmxlcyB1cCB0byB6IGJpdHNcbiAgICAgICAgICAgICAgICBpZigoZiA8PD0gMSkgPD0gdGhpcy5jWysreHBdKVxuICAgICAgICAgICAgICAgICAgYnJlYWs7ICAgICAgICAgICAgICAvLyBlbm91Z2ggY29kZXMgdG8gdXNlIHVwIGogYml0c1xuICAgICAgICAgICAgICAgIGYgLT0gdGhpcy5jW3hwXTsgICAgICAgICAgIC8vIGVsc2UgZGVkdWN0IGNvZGVzIGZyb20gcGF0dGVybnNcbiAgICAgICAgICAgICAgfVxuXHQgICAgfVxuICAgICAgICAgIH1cbiAgICAgICAgICB6ID0gMSA8PCBqOyAgICAgICAgICAgICAgICAgLy8gdGFibGUgZW50cmllcyBmb3Igai1iaXQgdGFibGVcblxuXHQgIC8vIGFsbG9jYXRlIG5ldyB0YWJsZVxuICAgICAgICAgIGlmICh0aGlzLmhuWzBdICsgeiA+IE1BTlkpeyAgICAgICAvLyAobm90ZTogZG9lc24ndCBtYXR0ZXIgZm9yIGZpeGVkKVxuICAgICAgICAgICAgcmV0dXJuIFpfREFUQV9FUlJPUjsgICAgICAgLy8gb3ZlcmZsb3cgb2YgTUFOWVxuICAgICAgICAgIH1cbiAgICAgICAgICB0aGlzLnVbaF0gPSBxID0gLypocCsqLyB0aGlzLmhuWzBdOyAgIC8vIERFQlVHXG4gICAgICAgICAgdGhpcy5oblswXSArPSB6O1xuIFxuXHQgIC8vIGNvbm5lY3QgdG8gbGFzdCB0YWJsZSwgaWYgdGhlcmUgaXMgb25lXG5cdCAgaWYoaCE9MCl7XG4gICAgICAgICAgICB0aGlzLnhbaF09aTsgICAgICAgICAgIC8vIHNhdmUgcGF0dGVybiBmb3IgYmFja2luZyB1cFxuICAgICAgICAgICAgdGhpcy5yWzBdPWo7ICAgICAvLyBiaXRzIGluIHRoaXMgdGFibGVcbiAgICAgICAgICAgIHRoaXMuclsxXT1sOyAgICAgLy8gYml0cyB0byBkdW1wIGJlZm9yZSB0aGlzIHRhYmxlXG4gICAgICAgICAgICBqPWk+Pj4odyAtIGwpO1xuICAgICAgICAgICAgdGhpcy5yWzJdID0gKHEgLSB0aGlzLnVbaC0xXSAtIGopOyAgICAgICAgICAgICAgIC8vIG9mZnNldCB0byB0aGlzIHRhYmxlXG4gICAgICAgICAgICBhcnJheUNvcHkodGhpcy5yLCAwLCBocCwgKHRoaXMudVtoLTFdK2opKjMsIDMpOyAvLyBjb25uZWN0IHRvIGxhc3QgdGFibGVcbiAgICAgICAgICB9XG4gICAgICAgICAgZWxzZXtcbiAgICAgICAgICAgIHRbMF0gPSBxOyAgICAgICAgICAgICAgIC8vIGZpcnN0IHRhYmxlIGlzIHJldHVybmVkIHJlc3VsdFxuXHQgIH1cbiAgICAgICAgfVxuXG5cdC8vIHNldCB1cCB0YWJsZSBlbnRyeSBpbiByXG4gICAgICAgIHRoaXMuclsxXSA9IChrIC0gdyk7XG4gICAgICAgIGlmIChwID49IG4pe1xuICAgICAgICAgIHRoaXMuclswXSA9IDEyOCArIDY0OyAgICAgIC8vIG91dCBvZiB2YWx1ZXMtLWludmFsaWQgY29kZVxuXHR9XG4gICAgICAgIGVsc2UgaWYgKHZbcF0gPCBzKXtcbiAgICAgICAgICB0aGlzLnJbMF0gPSAodGhpcy52W3BdIDwgMjU2ID8gMCA6IDMyICsgNjQpOyAgLy8gMjU2IGlzIGVuZC1vZi1ibG9ja1xuICAgICAgICAgIHRoaXMuclsyXSA9IHRoaXMudltwKytdOyAgICAgICAgICAvLyBzaW1wbGUgY29kZSBpcyBqdXN0IHRoZSB2YWx1ZVxuICAgICAgICB9XG4gICAgICAgIGVsc2V7XG4gICAgICAgICAgdGhpcy5yWzBdPShlW3RoaXMudltwXS1zXSsxNis2NCk7IC8vIG5vbi1zaW1wbGUtLWxvb2sgdXAgaW4gbGlzdHNcbiAgICAgICAgICB0aGlzLnJbMl09ZFt0aGlzLnZbcCsrXSAtIHNdO1xuICAgICAgICB9XG5cbiAgICAgICAgLy8gZmlsbCBjb2RlLWxpa2UgZW50cmllcyB3aXRoIHJcbiAgICAgICAgZj0xPDwoay13KTtcbiAgICAgICAgZm9yIChqPWk+Pj53O2o8ejtqKz1mKXtcbiAgICAgICAgICBhcnJheUNvcHkodGhpcy5yLCAwLCBocCwgKHEraikqMywgMyk7XG5cdH1cblxuXHQvLyBiYWNrd2FyZHMgaW5jcmVtZW50IHRoZSBrLWJpdCBjb2RlIGlcbiAgICAgICAgZm9yIChqID0gMSA8PCAoayAtIDEpOyAoaSAmIGopIT0wOyBqID4+Pj0gMSl7XG4gICAgICAgICAgaSBePSBqO1xuXHR9XG4gICAgICAgIGkgXj0gajtcblxuXHQvLyBiYWNrdXAgb3ZlciBmaW5pc2hlZCB0YWJsZXNcbiAgICAgICAgbWFzayA9ICgxIDw8IHcpIC0gMTsgICAgICAvLyBuZWVkZWQgb24gSFAsIGNjIC1PIGJ1Z1xuICAgICAgICB3aGlsZSAoKGkgJiBtYXNrKSAhPSB0aGlzLnhbaF0pe1xuICAgICAgICAgIGgtLTsgICAgICAgICAgICAgICAgICAgIC8vIGRvbid0IG5lZWQgdG8gdXBkYXRlIHFcbiAgICAgICAgICB3IC09IGw7XG4gICAgICAgICAgbWFzayA9ICgxIDw8IHcpIC0gMTtcbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgICAvLyBSZXR1cm4gWl9CVUZfRVJST1IgaWYgd2Ugd2VyZSBnaXZlbiBhbiBpbmNvbXBsZXRlIHRhYmxlXG4gICAgcmV0dXJuIHkgIT0gMCAmJiBnICE9IDEgPyBaX0JVRl9FUlJPUiA6IFpfT0s7XG59XG5cbkluZlRyZWUucHJvdG90eXBlLmluZmxhdGVfdHJlZXNfYml0cyA9IGZ1bmN0aW9uKGMsIGJiLCB0YiwgaHAsIHopIHtcbiAgICB2YXIgcmVzdWx0O1xuICAgIHRoaXMuaW5pdFdvcmtBcmVhKDE5KTtcbiAgICB0aGlzLmhuWzBdPTA7XG4gICAgcmVzdWx0ID0gdGhpcy5odWZ0X2J1aWxkKGMsIDAsIDE5LCAxOSwgbnVsbCwgbnVsbCwgdGIsIGJiLCBocCwgdGhpcy5obiwgdGhpcy52KTtcblxuICAgIGlmKHJlc3VsdCA9PSBaX0RBVEFfRVJST1Ipe1xuICAgICAgei5tc2cgPSBcIm92ZXJzdWJzY3JpYmVkIGR5bmFtaWMgYml0IGxlbmd0aHMgdHJlZVwiO1xuICAgIH1cbiAgICBlbHNlIGlmKHJlc3VsdCA9PSBaX0JVRl9FUlJPUiB8fCBiYlswXSA9PSAwKXtcbiAgICAgIHoubXNnID0gXCJpbmNvbXBsZXRlIGR5bmFtaWMgYml0IGxlbmd0aHMgdHJlZVwiO1xuICAgICAgcmVzdWx0ID0gWl9EQVRBX0VSUk9SO1xuICAgIH1cbiAgICByZXR1cm4gcmVzdWx0O1xufVxuXG5JbmZUcmVlLnByb3RvdHlwZS5pbmZsYXRlX3RyZWVzX2R5bmFtaWMgPSBmdW5jdGlvbihubCwgbmQsIGMsIGJsLCBiZCwgdGwsIHRkLCBocCwgeikge1xuICAgIHZhciByZXN1bHQ7XG5cbiAgICAvLyBidWlsZCBsaXRlcmFsL2xlbmd0aCB0cmVlXG4gICAgdGhpcy5pbml0V29ya0FyZWEoMjg4KTtcbiAgICB0aGlzLmhuWzBdPTA7XG4gICAgcmVzdWx0ID0gdGhpcy5odWZ0X2J1aWxkKGMsIDAsIG5sLCAyNTcsIGNwbGVucywgY3BsZXh0LCB0bCwgYmwsIGhwLCB0aGlzLmhuLCB0aGlzLnYpO1xuICAgIGlmIChyZXN1bHQgIT0gWl9PSyB8fCBibFswXSA9PSAwKXtcbiAgICAgIGlmKHJlc3VsdCA9PSBaX0RBVEFfRVJST1Ipe1xuICAgICAgICB6Lm1zZyA9IFwib3ZlcnN1YnNjcmliZWQgbGl0ZXJhbC9sZW5ndGggdHJlZVwiO1xuICAgICAgfVxuICAgICAgZWxzZSBpZiAocmVzdWx0ICE9IFpfTUVNX0VSUk9SKXtcbiAgICAgICAgei5tc2cgPSBcImluY29tcGxldGUgbGl0ZXJhbC9sZW5ndGggdHJlZVwiO1xuICAgICAgICByZXN1bHQgPSBaX0RBVEFfRVJST1I7XG4gICAgICB9XG4gICAgICByZXR1cm4gcmVzdWx0O1xuICAgIH1cblxuICAgIC8vIGJ1aWxkIGRpc3RhbmNlIHRyZWVcbiAgICB0aGlzLmluaXRXb3JrQXJlYSgyODgpO1xuICAgIHJlc3VsdCA9IHRoaXMuaHVmdF9idWlsZChjLCBubCwgbmQsIDAsIGNwZGlzdCwgY3BkZXh0LCB0ZCwgYmQsIGhwLCB0aGlzLmhuLCB0aGlzLnYpO1xuXG4gICAgaWYgKHJlc3VsdCAhPSBaX09LIHx8IChiZFswXSA9PSAwICYmIG5sID4gMjU3KSl7XG4gICAgICBpZiAocmVzdWx0ID09IFpfREFUQV9FUlJPUil7XG4gICAgICAgIHoubXNnID0gXCJvdmVyc3Vic2NyaWJlZCBkaXN0YW5jZSB0cmVlXCI7XG4gICAgICB9XG4gICAgICBlbHNlIGlmIChyZXN1bHQgPT0gWl9CVUZfRVJST1IpIHtcbiAgICAgICAgei5tc2cgPSBcImluY29tcGxldGUgZGlzdGFuY2UgdHJlZVwiO1xuICAgICAgICByZXN1bHQgPSBaX0RBVEFfRVJST1I7XG4gICAgICB9XG4gICAgICBlbHNlIGlmIChyZXN1bHQgIT0gWl9NRU1fRVJST1Ipe1xuICAgICAgICB6Lm1zZyA9IFwiZW1wdHkgZGlzdGFuY2UgdHJlZSB3aXRoIGxlbmd0aHNcIjtcbiAgICAgICAgcmVzdWx0ID0gWl9EQVRBX0VSUk9SO1xuICAgICAgfVxuICAgICAgcmV0dXJuIHJlc3VsdDtcbiAgICB9XG5cbiAgICByZXR1cm4gWl9PSztcbn1cbi8qXG4gIHN0YXRpYyBpbnQgaW5mbGF0ZV90cmVlc19maXhlZChpbnRbXSBibCwgIC8vbGl0ZXJhbCBkZXNpcmVkL2FjdHVhbCBiaXQgZGVwdGhcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGludFtdIGJkLCAgLy9kaXN0YW5jZSBkZXNpcmVkL2FjdHVhbCBiaXQgZGVwdGhcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGludFtdW10gdGwsLy9saXRlcmFsL2xlbmd0aCB0cmVlIHJlc3VsdFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW50W11bXSB0ZCwvL2Rpc3RhbmNlIHRyZWUgcmVzdWx0IFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgWlN0cmVhbSB6ICAvL2ZvciBtZW1vcnkgYWxsb2NhdGlvblxuXHRcdFx0XHQgKXtcblxuKi9cblxuZnVuY3Rpb24gaW5mbGF0ZV90cmVlc19maXhlZChibCwgYmQsIHRsLCB0ZCwgeikge1xuICAgIGJsWzBdPWZpeGVkX2JsO1xuICAgIGJkWzBdPWZpeGVkX2JkO1xuICAgIHRsWzBdPWZpeGVkX3RsO1xuICAgIHRkWzBdPWZpeGVkX3RkO1xuICAgIHJldHVybiBaX09LO1xufVxuXG5JbmZUcmVlLnByb3RvdHlwZS5pbml0V29ya0FyZWEgPSBmdW5jdGlvbih2c2l6ZSl7XG4gICAgaWYodGhpcy5obj09bnVsbCl7XG4gICAgICAgIHRoaXMuaG49bmV3IEludDMyQXJyYXkoMSk7XG4gICAgICAgIHRoaXMudj1uZXcgSW50MzJBcnJheSh2c2l6ZSk7XG4gICAgICAgIHRoaXMuYz1uZXcgSW50MzJBcnJheShCTUFYKzEpO1xuICAgICAgICB0aGlzLnI9bmV3IEludDMyQXJyYXkoMyk7XG4gICAgICAgIHRoaXMudT1uZXcgSW50MzJBcnJheShCTUFYKTtcbiAgICAgICAgdGhpcy54PW5ldyBJbnQzMkFycmF5KEJNQVgrMSk7XG4gICAgfVxuICAgIGlmKHRoaXMudi5sZW5ndGg8dnNpemUpeyBcbiAgICAgICAgdGhpcy52PW5ldyBJbnQzMkFycmF5KHZzaXplKTsgXG4gICAgfVxuICAgIGZvcih2YXIgaT0wOyBpPHZzaXplOyBpKyspe3RoaXMudltpXT0wO31cbiAgICBmb3IodmFyIGk9MDsgaTxCTUFYKzE7IGkrKyl7dGhpcy5jW2ldPTA7fVxuICAgIGZvcih2YXIgaT0wOyBpPDM7IGkrKyl7dGhpcy5yW2ldPTA7fVxuLy8gIGZvcihpbnQgaT0wOyBpPEJNQVg7IGkrKyl7dVtpXT0wO31cbiAgICBhcnJheUNvcHkodGhpcy5jLCAwLCB0aGlzLnUsIDAsIEJNQVgpO1xuLy8gIGZvcihpbnQgaT0wOyBpPEJNQVgrMTsgaSsrKXt4W2ldPTA7fVxuICAgIGFycmF5Q29weSh0aGlzLmMsIDAsIHRoaXMueCwgMCwgQk1BWCsxKTtcbn1cblxudmFyIHRlc3RBcnJheSA9IG5ldyBVaW50OEFycmF5KDEpO1xudmFyIGhhc1N1YmFycmF5ID0gKHR5cGVvZiB0ZXN0QXJyYXkuc3ViYXJyYXkgPT09ICdmdW5jdGlvbicpO1xudmFyIGhhc1NsaWNlID0gZmFsc2U7IC8qICh0eXBlb2YgdGVzdEFycmF5LnNsaWNlID09PSAnZnVuY3Rpb24nKTsgKi8gLy8gQ2hyb21lIHNsaWNlIHBlcmZvcm1hbmNlIGlzIHNvIGRpcmUgdGhhdCB3ZSdyZSBjdXJyZW50bHkgbm90IHVzaW5nIGl0Li4uXG5cbmZ1bmN0aW9uIGFycmF5Q29weShzcmMsIHNyY09mZnNldCwgZGVzdCwgZGVzdE9mZnNldCwgY291bnQpIHtcbiAgICBpZiAoY291bnQgPT0gMCkge1xuICAgICAgICByZXR1cm47XG4gICAgfSBcbiAgICBpZiAoIXNyYykge1xuICAgICAgICB0aHJvdyBcIlVuZGVmIHNyY1wiO1xuICAgIH0gZWxzZSBpZiAoIWRlc3QpIHtcbiAgICAgICAgdGhyb3cgXCJVbmRlZiBkZXN0XCI7XG4gICAgfVxuXG4gICAgaWYgKHNyY09mZnNldCA9PSAwICYmIGNvdW50ID09IHNyYy5sZW5ndGgpIHtcbiAgICAgICAgYXJyYXlDb3B5X2Zhc3Qoc3JjLCBkZXN0LCBkZXN0T2Zmc2V0KTtcbiAgICB9IGVsc2UgaWYgKGhhc1N1YmFycmF5KSB7XG4gICAgICAgIGFycmF5Q29weV9mYXN0KHNyYy5zdWJhcnJheShzcmNPZmZzZXQsIHNyY09mZnNldCArIGNvdW50KSwgZGVzdCwgZGVzdE9mZnNldCk7IFxuICAgIH0gZWxzZSBpZiAoc3JjLkJZVEVTX1BFUl9FTEVNRU5UID09IDEgJiYgY291bnQgPiAxMDApIHtcbiAgICAgICAgYXJyYXlDb3B5X2Zhc3QobmV3IFVpbnQ4QXJyYXkoc3JjLmJ1ZmZlciwgc3JjLmJ5dGVPZmZzZXQgKyBzcmNPZmZzZXQsIGNvdW50KSwgZGVzdCwgZGVzdE9mZnNldCk7XG4gICAgfSBlbHNlIHsgXG4gICAgICAgIGFycmF5Q29weV9zbG93KHNyYywgc3JjT2Zmc2V0LCBkZXN0LCBkZXN0T2Zmc2V0LCBjb3VudCk7XG4gICAgfVxuXG59XG5cbmZ1bmN0aW9uIGFycmF5Q29weV9zbG93KHNyYywgc3JjT2Zmc2V0LCBkZXN0LCBkZXN0T2Zmc2V0LCBjb3VudCkge1xuXG4gICAgLy8gZGxvZygnX3Nsb3cgY2FsbDogc3JjT2Zmc2V0PScgKyBzcmNPZmZzZXQgKyAnOyBkZXN0T2Zmc2V0PScgKyBkZXN0T2Zmc2V0ICsgJzsgY291bnQ9JyArIGNvdW50KTtcblxuICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGNvdW50OyArK2kpIHtcbiAgICAgICAgZGVzdFtkZXN0T2Zmc2V0ICsgaV0gPSBzcmNbc3JjT2Zmc2V0ICsgaV07XG4gICAgfVxufVxuXG5mdW5jdGlvbiBhcnJheUNvcHlfZmFzdChzcmMsIGRlc3QsIGRlc3RPZmZzZXQpIHtcbiAgICBkZXN0LnNldChzcmMsIGRlc3RPZmZzZXQpO1xufVxuXG5cbiAgLy8gbGFyZ2VzdCBwcmltZSBzbWFsbGVyIHRoYW4gNjU1MzZcbnZhciBBRExFUl9CQVNFPTY1NTIxOyBcbiAgLy8gTk1BWCBpcyB0aGUgbGFyZ2VzdCBuIHN1Y2ggdGhhdCAyNTVuKG4rMSkvMiArIChuKzEpKEJBU0UtMSkgPD0gMl4zMi0xXG52YXIgQURMRVJfTk1BWD01NTUyO1xuXG5mdW5jdGlvbiBhZGxlcjMyKGFkbGVyLCAvKiBieXRlW10gKi8gYnVmLCAgaW5kZXgsIGxlbil7XG4gICAgaWYoYnVmID09IG51bGwpeyByZXR1cm4gMTsgfVxuXG4gICAgdmFyIHMxPWFkbGVyJjB4ZmZmZjtcbiAgICB2YXIgczI9KGFkbGVyPj4xNikmMHhmZmZmO1xuICAgIHZhciBrO1xuXG4gICAgd2hpbGUobGVuID4gMCkge1xuICAgICAgaz1sZW48QURMRVJfTk1BWD9sZW46QURMRVJfTk1BWDtcbiAgICAgIGxlbi09aztcbiAgICAgIHdoaWxlKGs+PTE2KXtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIGstPTE2O1xuICAgICAgfVxuICAgICAgaWYoayE9MCl7XG4gICAgICAgIGRve1xuICAgICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICB9XG4gICAgICAgIHdoaWxlKC0tayE9MCk7XG4gICAgICB9XG4gICAgICBzMSU9QURMRVJfQkFTRTtcbiAgICAgIHMyJT1BRExFUl9CQVNFO1xuICAgIH1cbiAgICByZXR1cm4gKHMyPDwxNil8czE7XG59XG5cblxuXG5mdW5jdGlvbiBqc3psaWJfaW5mbGF0ZV9idWZmZXIoYnVmZmVyLCBzdGFydCwgbGVuZ3RoLCBhZnRlclVuY09mZnNldCkge1xuICAgIGlmICghc3RhcnQpIHtcbiAgICAgICAgYnVmZmVyID0gbmV3IFVpbnQ4QXJyYXkoYnVmZmVyKTtcbiAgICB9IGVsc2UgaWYgKCFsZW5ndGgpIHtcbiAgICAgICAgYnVmZmVyID0gbmV3IFVpbnQ4QXJyYXkoYnVmZmVyLCBzdGFydCwgYnVmZmVyLmJ5dGVMZW5ndGggLSBzdGFydCk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgYnVmZmVyID0gbmV3IFVpbnQ4QXJyYXkoYnVmZmVyLCBzdGFydCwgbGVuZ3RoKTtcbiAgICB9XG5cbiAgICB2YXIgeiA9IG5ldyBaU3RyZWFtKCk7XG4gICAgei5pbmZsYXRlSW5pdChERUZfV0JJVFMsIHRydWUpO1xuICAgIHoubmV4dF9pbiA9IGJ1ZmZlcjtcbiAgICB6Lm5leHRfaW5faW5kZXggPSAwO1xuICAgIHouYXZhaWxfaW4gPSBidWZmZXIubGVuZ3RoO1xuXG4gICAgdmFyIG9CbG9ja0xpc3QgPSBbXTtcbiAgICB2YXIgdG90YWxTaXplID0gMDtcbiAgICB3aGlsZSAodHJ1ZSkge1xuICAgICAgICB2YXIgb2J1ZiA9IG5ldyBVaW50OEFycmF5KDMyMDAwKTtcbiAgICAgICAgei5uZXh0X291dCA9IG9idWY7XG4gICAgICAgIHoubmV4dF9vdXRfaW5kZXggPSAwO1xuICAgICAgICB6LmF2YWlsX291dCA9IG9idWYubGVuZ3RoO1xuICAgICAgICB2YXIgc3RhdHVzID0gei5pbmZsYXRlKFpfTk9fRkxVU0gpO1xuICAgICAgICBpZiAoc3RhdHVzICE9IFpfT0sgJiYgc3RhdHVzICE9IFpfU1RSRUFNX0VORCAmJiBzdGF0dXMgIT0gWl9CVUZfRVJST1IpIHtcbiAgICAgICAgICAgIHRocm93IHoubXNnO1xuICAgICAgICB9XG4gICAgICAgIGlmICh6LmF2YWlsX291dCAhPSAwKSB7XG4gICAgICAgICAgICB2YXIgbmV3b2IgPSBuZXcgVWludDhBcnJheShvYnVmLmxlbmd0aCAtIHouYXZhaWxfb3V0KTtcbiAgICAgICAgICAgIGFycmF5Q29weShvYnVmLCAwLCBuZXdvYiwgMCwgKG9idWYubGVuZ3RoIC0gei5hdmFpbF9vdXQpKTtcbiAgICAgICAgICAgIG9idWYgPSBuZXdvYjtcbiAgICAgICAgfVxuICAgICAgICBvQmxvY2tMaXN0LnB1c2gob2J1Zik7XG4gICAgICAgIHRvdGFsU2l6ZSArPSBvYnVmLmxlbmd0aDtcbiAgICAgICAgaWYgKHN0YXR1cyA9PSBaX1NUUkVBTV9FTkQgfHwgc3RhdHVzID09IFpfQlVGX0VSUk9SKSB7XG4gICAgICAgICAgICBicmVhaztcbiAgICAgICAgfVxuICAgIH1cblxuICAgIGlmIChhZnRlclVuY09mZnNldCkge1xuICAgICAgICBhZnRlclVuY09mZnNldFswXSA9IChzdGFydCB8fCAwKSArIHoubmV4dF9pbl9pbmRleDtcbiAgICB9XG5cbiAgICBpZiAob0Jsb2NrTGlzdC5sZW5ndGggPT0gMSkge1xuICAgICAgICByZXR1cm4gb0Jsb2NrTGlzdFswXS5idWZmZXI7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdmFyIG91dCA9IG5ldyBVaW50OEFycmF5KHRvdGFsU2l6ZSk7XG4gICAgICAgIHZhciBjdXJzb3IgPSAwO1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG9CbG9ja0xpc3QubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgIHZhciBiID0gb0Jsb2NrTGlzdFtpXTtcbiAgICAgICAgICAgIGFycmF5Q29weShiLCAwLCBvdXQsIGN1cnNvciwgYi5sZW5ndGgpO1xuICAgICAgICAgICAgY3Vyc29yICs9IGIubGVuZ3RoO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBvdXQuYnVmZmVyO1xuICAgIH1cbn1cblxuaWYgKHR5cGVvZihtb2R1bGUpICE9PSAndW5kZWZpbmVkJykge1xuICBtb2R1bGUuZXhwb3J0cyA9IHtcbiAgICBpbmZsYXRlQnVmZmVyOiBqc3psaWJfaW5mbGF0ZV9idWZmZXIsXG4gICAgYXJyYXlDb3B5OiBhcnJheUNvcHlcbiAgfTtcbn1cbiJdfQ==
