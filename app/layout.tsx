import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "RDP Web",
  description:
    "A local-first, WebAssembly recombination detection and review workbench for aligned nucleotide sequences.",
  other: {
    "codex-preview": "development",
  },
  icons: {
    icon: "./favicon.svg",
    shortcut: "./favicon.svg",
  },
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
