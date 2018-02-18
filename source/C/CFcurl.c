#include "CFcurl.h"

static size_t write_data(void *ptr, size_t size, size_t nmemb, void *stream) {
  size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
  return written;
  }

void OXCALL fGet OXARGS {
    OxLibCheckType(OX_STRING,pv,0,1);
    OxLibCheckType(OX_INT,pv,0,2);
    CURL *curl_handle;
    FILE *pagefile;
    CURLcode res;

    curl_global_init(CURL_GLOBAL_ALL);
    curl_handle = curl_easy_init();
    curl_easy_setopt(curl_handle, CURLOPT_URL, OxStr(pv,0));

    curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 1L);
    curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);

    OxSetInt(rtn,0,FALSE);
    pagefile = fopen(OxStr(pv,1), "wb");
    if(pagefile) {
        curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, pagefile);
        res = curl_easy_perform(curl_handle);
        if(res != CURLE_OK) {
            if (OxInt(pv,2))
                fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));
            }
        else
            OxSetInt(rtn,0,TRUE);
        fclose(pagefile);
        }
    else {
        if (OxInt(pv,2)) fprintf(stderr, "file creation failed: %s\n");
        }
    curl_easy_cleanup(curl_handle);
  }

void OXCALL fPut OXARGS {
  OxLibCheckType(OX_STRING,pv,0,1);
  OxLibCheckType(OX_INT,pv,0,2);
  CURL *curl;
  CURLcode res;
  struct stat file_info;
  double speed_upload, total_time;
  FILE *fd;

  fd = fopen(OxStr(pv,1), "rb"); /* open file to upload */
  OxSetInt(rtn,0,FALSE); /* fail unless we get to end */
  if (!fd) return;
  if (fstat(fileno(fd), &file_info) != 0) return;
  curl = curl_easy_init();
  if (curl) {
    curl_easy_setopt(curl, CURLOPT_URL,OxStr(pv,0));
    curl_easy_setopt(curl, CURLOPT_UPLOAD, 1L);
    curl_easy_setopt(curl, CURLOPT_READDATA, fd);
    curl_easy_setopt(curl, CURLOPT_INFILESIZE_LARGE,
                     (curl_off_t)file_info.st_size);
    curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
    res = curl_easy_perform(curl);
    if(res != CURLE_OK) {
      if (OxInt(pv,2)) fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));
      }
    else { /* now extract transfer info */
      curl_easy_getinfo(curl, CURLINFO_SPEED_UPLOAD, &speed_upload);
      curl_easy_getinfo(curl, CURLINFO_TOTAL_TIME, &total_time);
      OxSetInt(rtn,0,TRUE);
      if (OxInt(pv,2))
        fprintf(stderr, "Speed: %.3f bytes/sec during %.3f seconds\n",speed_upload, total_time);
    }
  curl_easy_cleanup(curl);
  }
  fclose(fd);
  }

void OXCALL fPost OXARGS {
  OxLibCheckType(OX_STRING,pv,0,1);
  OxLibCheckType(OX_INT,pv,0,2);
  CURL *curl;
  CURLcode res;
  curl_global_init(CURL_GLOBAL_ALL);
  curl = curl_easy_init();
  OxSetInt(rtn,0,FALSE);
  if(curl) {
    /* set the URL that is about to receive our POST.  */
    curl_easy_setopt(curl, CURLOPT_URL, OxStr(pv,0) );
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, OxStr(pv,1) );
    res = curl_easy_perform(curl);
    /* Check for errors */
    if(res != CURLE_OK) {
      if (OxInt(pv,2))
            fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
        }
    else
        OxSetInt(rtn,0,TRUE);
    curl_easy_cleanup(curl);
  }
  curl_global_cleanup();
}
