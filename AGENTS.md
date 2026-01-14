# Repository Guidelines

## Project Structure

This repository is a GitHub Pages/Jekyll site for HWRS 504 (Fall 2025).

- `_lectures/`: lecture pages (e.g., `lecture-01.md`) rendered with Jekyll layouts.
- `_schedules/`: weekly schedule entries (e.g., `week-01.md`).
- `_announcements/`: announcement posts (e.g., `week-0.md`).
- `_layouts/`, `_includes/`, `_sass/`: site templates, partials, and styling.
- `assets/`: static files (images, calendars, PDFs, etc.).
- Top-level `*.md` (e.g., `schedule.md`, `policy.md`): site pages.
- `pluto_notebooks/`: Pluto/Julia notebooks and related figures (not required for site build).

## Build, Test, and Development Commands

Prereqs: Ruby + Bundler.

- `bundle install`: install Jekyll/GitHub Pages dependencies.
- `bundle exec jekyll serve`: run the site locally (default at `http://127.0.0.1:4000/`).
- `bundle exec jekyll serve --livereload`: rebuild on change with live reload.
- `bundle exec jekyll build`: produce the static site output in `_site/`.

## Coding Style & Naming Conventions

- Markdown pages should include YAML front matter (`---` block) when rendered by Jekyll.
- Keep filenames consistent with existing patterns: `lecture-XX.md`, `week-XX.md`.
- Use 2-space indentation in YAML and Liquid templates; keep lines readable and avoid trailing whitespace.
- Prefer relative links for in-repo assets (e.g., `assets/images/...`) and keep paths stable.

## Testing Guidelines

No automated test suite is configured. Validate changes by:

- Running `bundle exec jekyll build` to catch Liquid/front-matter errors.
- Spot-checking affected pages locally with `bundle exec jekyll serve`.

## Commit & Pull Request Guidelines

- Commits in this repo are typically short and action-oriented (e.g., “Updated module 17”, “Uploaded notebooks”).
- Write commit subjects in imperative mood: `Update lecture 12 slides`, `Add week-05 schedule`.
- PRs should describe what changed and why, link any related issue, and include screenshots for layout/styling changes.
- Avoid committing generated output (`_site/`) unless explicitly required for a special deployment workflow.

## Agent-Specific Notes

When making automated edits, keep diffs minimal, preserve existing page URLs/filenames, and verify the site builds locally.
