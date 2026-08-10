import type { NextConfig } from "next";

const githubPages = process.env.RDP_GITHUB_PAGES === "1";
const requestedBasePath = (process.env.RDP_BASE_PATH ?? "").trim();
const basePath = requestedBasePath
  ? `/${requestedBasePath.replace(/^\/+|\/+$/g, "")}`
  : "";

const nextConfig: NextConfig = {
  ...(githubPages
    ? {
        output: "export" as const,
        basePath,
        assetPrefix: basePath || undefined,
        trailingSlash: true,
        images: { unoptimized: true },
      }
    : {}),
};

export default nextConfig;
