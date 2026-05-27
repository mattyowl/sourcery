# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Sourcery is a web application for browsing and managing astronomical source catalogs. It serves a MongoDB-backed catalog as an interactive website with image overlays from various sky surveys, NED cross-matching, user-editable fields, and catalog export.

## Commands

**Install (development):**
```
python setup.py install --prefix=$HOME/local
```

**Build the MongoDB database from a FITS catalog:**
```
sourcery_build_db myconfig.yml
```

**Build/rebuild the image cache (downloads survey images per source):**
```
sourcery_build_cache myconfig.yml
```

**Run a local CherryPy development server (port 8080):**
```
sourcery_test myconfig.yml
# then open http://localhost:8080/sourcery
```

**Generate a password hash for user accounts:**
```
sourcery_password_hash
```

There is no automated test suite.

## Architecture

Everything is driven by a YAML configuration file (see `examples/E-D56/E-D56Clusters.yml` for the full annotated example). The config specifies the input FITS catalog, MongoDB collection name, image sources, editable fields, access control, and display options.

### Core data flow

1. `sourcery_build_db config.yml` — creates a `SourceBrowser` with `buildDatabase=True`, which reads the FITS catalog and writes it into MongoDB, optionally cross-matching against other FITS tables (`crossMatchCatalogs`).
2. `sourcery_build_cache config.yml` — calls `SourceBrowser.preprocess()`, which downloads/generates JPEG thumbnails for each source from web services (SDSS, DECaLS, PS1, unWISE) and local image directories (`imageDirs`, `tileDirs`), storing everything in `cacheDir`.
3. `sourcery_test config.yml` — starts CherryPy with `SourceBrowser` mounted at `/sourcery`. The web app serves the table and per-source pages directly from MongoDB and the image cache.

### Key classes and modules

- **`sourcery/sourceBrowser.py`** — The main class `SourceBrowser` is both a CherryPy application (all URL handlers are methods decorated with `@cherrypy.expose`) and the orchestrator for database building and cache generation. All web endpoints live here: `index()` (table view), `displaySourcePage()`, `downloadCatalog()`, `login()`/`logout()`, `updateMongoDB()`, etc.

- **`sourcery/catalogTools.py`** — Utility functions for image clipping, tan-resampling, coordinate matching. Standalone helpers that `sourceBrowser.py` imports.

- **`sourcery/tileDir.py`** — `TileDir` class handles surveys distributed as tile images (DES, KiDS, S82). Given a directory of `.jpg` tiles and a `WCSTab.fits` mapping tilenames to WCS, it stitches cutouts for each source.

- **`sourcery/sourceryAuth.py`** — CherryPy authentication tool (`check_auth`). Sessions use `pbkdf2_sha256` hashes. User roles are `editor` (can edit fields/classifications) or `viewer` (read-only). Auth is optional — if `userListFile` is absent in the config, the site is open.

- **`sourcery/specFeatures.py`** — Static list of spectral lines used by `makeSpectrumPlot`.

- **`bin/`** — Thin CLI wrappers that instantiate `SourceBrowser` with the appropriate flags and call the relevant method.

- **`cgi-bin/example-sourcery.wsgi`** — Entry point for Apache/mod_wsgi deployment (see `APACHE_DEPLOYMENT.md`).

### MongoDB structure

Two collections per deployment:
- `source` — the main catalog data (from the FITS file + cross-matches)
- A separate MongoDB collection (named by `MongoDBName` in the config) stores user-editable fields and classifications, cross-matched at runtime by position (`MongoDBCrossMatchRadiusArcmin`) so the input catalog can be swapped without losing user annotations.

### Image cache

All generated/fetched images go in `cacheDir` as JPEGs named by RA/dec (`%.5f_%.5f`). Survey images come from web APIs (SDSS, Legacy Survey/DECaLS, PS1, unWISE) or are cut from local FITS `imageDirs`. The cache is built once and incrementally updated for new sources only.

### Config-driven behaviour

Almost all features are opt-in via the YAML file: NED cross-matching (`addNEDMatches`), SDSS redshifts (`addSDSSRedshifts`), image sources (`addSDSSImage`, `addDECaLSImage`, …), contour overlays, user-editable fields, classifications, quick-link queries, and access control. The annotated example in `examples/E-D56/E-D56Clusters.yml` is the primary reference for available options.
