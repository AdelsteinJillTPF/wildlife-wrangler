/*
Use the occurrence data to evaluate the GAP range map for a species
*/

/* Load the GAP range csv, filter out some columns, rename others */
.mode csv
.import /users/nmtarr/Documents/ranges/indata/bYBCUx_CONUS_Range_2001v1.csv sp_range
.headers on
ALTER TABLE sp_range RENAME TO garb;
CREATE TABLE sp_range AS
                      SELECT strHUC12RNG AS HUC12RNG,
                             intGapOrigin AS intGAPOrigin,
                             intGapPres AS intGAPPresence,
                             intGapRepro AS intGAPReproduction,
                             intGapSeas AS intGAPSeason,
                             Origin AS strGAPOrigin,
                             Presence AS strGAPPresence,
                             Reproduction AS strGAPReproduction,
                             Season AS strGAPSeason
                      FROM garb;
DROP TABLE garb;

/*
Add a column to store gbif evaluation result - meaning, does the huc
completely contain an occurrence circle?  Circles that overlap huc boundaries
are left out of this, but could possibly be added later
*/
ALTER TABLE sp_range ADD COLUMN eval_gbif1 INTEGER;

/*  Select hucs that contain an occurrence circle */
CREATE TABLE blue AS
              SELECT shucs.HUC12RNG, shucs.geom_102008
              FROM shucs, occs
              WHERE Contains(shucs.geom_102008, occs.circle_albers);

/*  Record in sp_range that gap and gbif agreed on species presence */
UPDATE sp_range
SET eval_gbif1 = 1
WHERE EXISTS (SELECT HUC12RNG FROM blue WHERE HUC12RNG = sp_range.HUC12RNG);

/*  Find hucs that contained gbif occurrences, but were not in gaprange and
insert them into sp_range as new records */
INSERT INTO sp_range (HUC12RNG, eval_gbif1)
SELECT blue.HUC12RNG, 0
FROM blue LEFT JOIN sp_range ON sp_range.HUC12RNG = blue.HUC12RNG
WHERE sp_range.HUC12RNG IS NULL;

/*  For new records, put zeros in GAP range attribute fields  */
UPDATE sp_range
SET intGAPOrigin = 0,
    intGAPPresence = 0,
    intGAPReproduction = 0,
    intGAPSeason = 0
WHERE eval_gbif1 = 0;






/*  Create a version of sp_range with geometry  */
CREATE TABLE sp_geom AS
              SELECT sp_range.*, shucs.geom_102008
              FROM sp_range LEFT JOIN shucs ON sp_range.HUC12RNG = shucs.HUC12RNG;
SELECT RecoverGeometryColumn('sp_geom', 'geom_102008', 102008, 'POLYGON');

/* Export maps */
SELECT ExportSHP('sp_geom', 'geom_102008', 'sp_geom1', 'utf-8');
