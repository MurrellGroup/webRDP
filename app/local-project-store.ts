import { RdpProject } from "./rdp-core";

export interface AutosaveRecord {
  id: "current";
  savedAt: number;
  project: RdpProject;
}

const DATABASE = "rdp-web-local-projects";
const STORE = "autosaves";

function database(): Promise<IDBDatabase> {
  return new Promise((resolve, reject) => {
    const request = indexedDB.open(DATABASE, 1);
    request.onupgradeneeded = () => {
      if (!request.result.objectStoreNames.contains(STORE)) request.result.createObjectStore(STORE, { keyPath: "id" });
    };
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error ?? new Error("Could not open the local autosave database."));
  });
}

export async function loadAutosave(): Promise<AutosaveRecord | null> {
  const db = await database();
  try {
    return await new Promise((resolve, reject) => {
      const request = db.transaction(STORE, "readonly").objectStore(STORE).get("current");
      request.onsuccess = () => resolve((request.result as AutosaveRecord | undefined) ?? null);
      request.onerror = () => reject(request.error ?? new Error("Could not read the local autosave."));
    });
  } finally {
    db.close();
  }
}

export async function saveAutosave(project: RdpProject): Promise<number> {
  const db = await database();
  const savedAt = Date.now();
  try {
    await new Promise<void>((resolve, reject) => {
      const request = db.transaction(STORE, "readwrite").objectStore(STORE).put({ id: "current", savedAt, project } satisfies AutosaveRecord);
      request.onsuccess = () => resolve();
      request.onerror = () => reject(request.error ?? new Error("Could not save the local project."));
    });
    return savedAt;
  } finally {
    db.close();
  }
}

export async function clearAutosave(): Promise<void> {
  const db = await database();
  try {
    await new Promise<void>((resolve, reject) => {
      const request = db.transaction(STORE, "readwrite").objectStore(STORE).delete("current");
      request.onsuccess = () => resolve();
      request.onerror = () => reject(request.error ?? new Error("Could not clear the local autosave."));
    });
  } finally {
    db.close();
  }
}
