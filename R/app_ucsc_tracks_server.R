ucsc_track_server = function(root_dir, base_host_path, base_url){
  CFG_DIR = file.path(root_dir, ".track_configs")
  dir.create(CFG_DIR, showWarnings = FALSE)

  function(input, output, session) {
    rv = reactiveValues(
      bwDF = refresh_bw(root_dir),
      cfgDF = NULL,
      trackTxt = NULL
    )

    observe({
      req(rv$cfgDF)
      tdf = rv$cfgDF
      req(input$trackSize)
      req(input$viewLimitsMin)
      req(input$windowFun)
      tracks = sapply(seq_len(nrow(tdf)), function(i){
        x = as.character(tdf$file)[i]
        col = as.character(tdf$color)[i]
        makeTrack(
          base_host_path = base_host_path,
          base_url = base_url,
          file = x,
          type = input$fileType,
          name = basename(x),
          description = basename(x),
          size = input$trackSize,
          view_min = input$viewLimitsMin,
          view_max = input$viewLimitsMax,
          color = col,
          visibility = input$visibility,
          autoscale = input$autoScale,
          alwaysZero = input$alwaysZero,
          showZero = input$showZero,
          showYline = input$showYref,
          valYline = input$numYref,
          smoothingWindow = input$smoothWin,
          windowFun = input$windowFun)

      })
      rv$trackTxt = tracks
    })

    ## bw selection table
    output$bwTable = DT::renderDataTable({
      DT::datatable(rv$bwDF, filter = "top", options = list(pageLength = 25))
    })
    output$colorBy = renderUI(expr = {
      choices = colnames(rv$bwDF)#[grepl("V", colnames(rv$bwDF))]
      # sel = which(grepl("directory", choices))

      selectInput(inputId = "colorBy", label = "Color By", choices = choices, selected = choices[1])
    })
    output$colorAssign = renderUI(expr = {
      req(input$colorBy)
      req(rv$bwDF)
      req(input$bwTable_rows_all)
      grps = unique(rv$bwDF[as.numeric(input$bwTable_rows_all),][[input$colorBy]])

      gen_color_picker_ui(grps = grps, brew_name = input$colorStyle, is_free_color = input$freeColor)
    })

    observeEvent(input$removeSelected, {
      rvis = sort(as.integer(input$bwTable_rows_all))

      rsel = sort(as.integer(input$bwTable_rows_selected))
      rvis = setdiff(rvis, rsel)
      if(length(rsel) < 1){
        showNotification("Empty selection.", type = "error")
      }else{
        rv$bwDF = rv$bwDF[rvis,]
        showNotification(paste("Removed", length(rsel), "entries."))
      }

    })
    observeEvent(input$limitToSelected, {
      print("limitToSelected")
      rsel = sort(as.integer(input$bwTable_rows_selected))
      nr = nrow(rv$bwDF)
      if(length(rsel) < 1){
        showNotification("Empty selection.", type = "error")
      }else{
        rv$bwDF = rv$bwDF[rsel,]
        showNotification(paste("Removed", nr - length(rsel), "entries."))
      }

    })
    observeEvent(input$limitToVisible, {
      print("limitToVisible")
      rvis = sort(as.integer(input$bwTable_rows_all))
      nr = nrow(rv$bwDF)
      if(length(rvis) < 1){
        showNotification("Would be empty.", type = "error")
      }else if(length(rvis) == nr){
        showNotification("No filters.", type = "error")
      }else{
        rv$bwDF = rv$bwDF[rvis,]
        showNotification(paste("Removed", nr - length(rvis), "entries."))
      }

    })

    output$cfgLink = renderUI({
      req(rv$trackTxt)
      tmpf = tempfile(pattern = "tracks_", tmpdir = CFG_DIR)
      write.table(rv$trackTxt, file = tmpf, quote = F, row.names = F, col.names = F)
      tmp_url = sub(base_host_path, base_url, tmpf)
      tags$a(href = tmp_url, "to File", target="_blank")
    })

    output$ucscLink = renderUI({
      req(rv$trackTxt)
      tmpf = tempfile(pattern = "tracks_", tmpdir = CFG_DIR)
      write.table(rv$trackTxt, file = tmpf, quote = F, row.names = F, col.names = F)
      tmp_url = sub(base_host_path, base_url, tmpf)

      ucsc_URL = "https://genome.ucsc.edu/cgi-bin/hgTracks?hgct_customText="
      tags$a(href = paste0(ucsc_URL, tmp_url), "to UCSC", target="_blank")
    })

    output$cfgText = renderText({
      paste(rv$trackTxt, collapse = "\n")
    })

    output$cfgTable = DT::renderDataTable({
      req(rv$cfgDF)
      uniq_colors = as.character(unique(rv$cfgDF$color))
      tab_dt = DT::datatable(rv$cfgDF,
                             options = list(pageLength = 50))
      formatStyle(table = tab_dt,
                  "color",  backgroundColor = styleEqual(
                    levels = uniq_colors,
                    values = sapply(strsplit(uniq_colors, ","), function(x){
                      x = as.integer(x)/255
                      rgb(x[1], x[2], x[3])
                    })
                  ))
    })

    observeEvent({
      input$tabsetPanel
      input$trackSize
      input$viewLimitsMin
      input$viewLimitsMax
    },{

      if(input$tabsetPanel == "Track Configuration"){
        allr = as.numeric(input$bwTable_rows_all)

        req(input$colorBy)
        grps = unique(rv$bwDF[allr,][[input$colorBy]])

        colorKeys = as.character(rv$bwDF[allr, input$colorBy])
        track_colors = vapply(colorKeys, function(g){
          g = gsub("/", "_", g)
          input[[paste0("color_", g)]]
        }, FUN.VALUE = "a string")
        track_colors = hex2ucsc_rgb(track_colors)
        files = as.character(rv$bwDF[allr, ]$files)
        n = length(files)
        rv$cfgDF = data.frame(file = files,
                              track_height = rep(input$trackSize, n),
                              view_min = rep(input$viewLimitsMin, n),
                              view_max = rep(input$viewLimitsMax, n),
                              color = track_colors)
      }
    })

    observeEvent(
      eventExpr = {
        input$fileType
        input$refreshFiles
      },
      handlerExpr = {
        print(input$fileType)
        switch(input$fileType,
               bigWig = {
                 rv$bwDF = refresh_bw(root_dir)
               },
               bigBed = {
                 rv$bwDF = refresh_bb(root_dir)
               },
               stop("file type not recognized"))
      }
    )

    output$dlTracks = downloadHandler(filename = "tracks.txt", content = function(file){
      write.table(rv$trackTxt, file = file, quote = F, row.names = F, col.names = F)
    })
  }
}
