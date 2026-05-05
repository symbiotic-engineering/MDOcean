FROM python:3.11-slim

ENV PLAYWRIGHT_BROWSERS_PATH=/ms-playwright

RUN pip install --no-cache-dir playwright \
 && python -m playwright install --with-deps chromium \
 && chmod -R 777 /ms-playwright

WORKDIR /work