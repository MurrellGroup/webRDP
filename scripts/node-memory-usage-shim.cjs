// Some restricted build sandboxes expose Node without /proc, which makes
// process.memoryUsage() throw before Next can compile. GitHub-hosted runners do
// not need this fallback, but preloading it is harmless there.
const original = process.memoryUsage.bind(process);

function zeroUsage() {
  return {
    rss: 0,
    heapTotal: 0,
    heapUsed: 0,
    external: 0,
    arrayBuffers: 0,
  };
}

function safeMemoryUsage() {
  try {
    return original();
  } catch {
    return zeroUsage();
  }
}

safeMemoryUsage.rss = () => {
  try {
    return typeof original.rss === "function" ? original.rss() : original().rss;
  } catch {
    return 0;
  }
};

process.memoryUsage = safeMemoryUsage;

