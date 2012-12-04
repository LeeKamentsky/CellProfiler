import re
import wx
import wx.grid
import os

from cellprofiler.gui.regexp_editor import edit_regexp


class LDFile(object):
    def __init__(self, path):
        self.path = path
        self.filename = os.path.split(path)[1]
        self.channel = None
        self.metadata = {}

class LDFileList(object):
    def __init__(self, files = None, order=None):
        self.files = files or {}
        self.order = order or list(sorted(self.files.keys()))
        
    def __len__(self):
        return len(self.order)
    
    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return [self[i] for i in range(idx.start, idx.stop, idx.step)
                    if idx < len(self)]
        if idx < 0:
            idx = len(self) + idx
        return self.files[self.order[idx]]
    
    def __delitem__(self, idx):
        if not hasattr(idx, "__iter__"):
            idx = [idx]
        to_check = set()
        for i in reversed(sorted(idx)):
            if i > len(self.order):
                continue
            to_check.add(self.order[i])
            del self.order[i]
        all_paths = set(self.order)
        for path in to_check:
            if path not in all_paths:
                del self.files[path]
            
    
    def get_any_file(self):
        return self.files.values()[0]
    
    def add_files(self, files):
        for f in files:
            if f.path not in self.files:
                self.files[f.path] = f
                self.order.append(f.path)
                
    def remove_files(self, paths):
        for path in paths:
            del self.files[path]
        self.order = filter((lambda x: x in self.files), order)
        
    def remove_selection(self, selection):
        self.remove_files(selection.get_selection())
        
    def sort(self):
        self.order = list(sorted(self.files.keys()))
            
    def get_files(self):
        return [self.files[p] for p in self.order]
    
    def get_selection(self, selection):
        return LDFileList(dict([(path, self.files[path]) 
                                for path in selection.get_selection()
                                if path in self.files]))
    
    def get_filtered_files(self, filter_fn):
        return LDSelection([path for path in self.files
                            if filter_fn(self.files[path])])
    
    def get_metadata_keys(self):
        keys = set()
        for f in self.files.values():
            keys.update(f.metadata.keys())
        return sorted(keys)
    
class LDSelection(object):
    def __init__(self, paths = []):
        self.paths = set(paths)
        
    def __len__(self):
        return len(self.paths)
    
    def __iter__(self):
        for path in self.paths:
            yield path
        
    def get_selection(self):
        return sorted(self.paths)
    
    def add_to_selection(self, paths):
        self.paths.add(paths)
        
    def remove_from_selection(self, paths):
        self.paths.difference_update(paths)

class LDDropTarget(wx.PyDropTarget):
    def __init__(self, file_callback_fn, text_callback_fn, 
                 over_callback_fn=None,
                 leave_callback_fn=None):
        super(self.__class__, self).__init__()
        self.file_callback_fn = file_callback_fn
        self.text_callback_fn = text_callback_fn
        self.over_callback_fn = over_callback_fn
        self.leave_callback_fn = leave_callback_fn
        self.file_data_object = wx.FileDataObject()
        self.text_data_object = wx.TextDataObject()
        self.composite_data_object = wx.DataObjectComposite()
        self.composite_data_object.Add(self.file_data_object, True)
        self.composite_data_object.Add(self.text_data_object)
        self.SetDataObject(self.composite_data_object)
        
    def OnDropFiles(self, x, y, filenames):
        self.file_callback_fn(x, y, filenames)
        
    def OnDropText(self, x, y, text):
        self.text_callback_fn(x, y, text)
        
    def OnEnter(self, x, y, d):
        return wx.DragCopy
    
    def OnLeave(self):
        if self.leave_callback_fn is not None:
            self.leave_callback_fn()
        return wx.PyDropTarget.OnLeave(self)
        
    def OnDragOver(self, x, y, d):
        if self.over_callback_fn is not None:
            return self.over_callback_fn(x, y, d)
        return wx.DragCopy
    
    def OnData(self, x, y, d):
        if self.GetData():
            df = self.composite_data_object.GetReceivedFormat().GetType()
            if  df in (wx.DF_TEXT, wx.DF_UNICODETEXT):
                self.OnDropText(x, y, self.text_data_object.GetText())
            elif df == wx.DF_FILENAME:
                self.OnDropFiles(x, y,
                                 self.file_data_object.GetFilenames())
        return wx.DragCopy
        
    def OnDrop(self, x, y):
        return True
    
class LDSelectionUIController(object):
    def __init__(self, panel, file_list_controller):
        self.file_list_controller = file_list_controller
        self.panel = panel
        assert isinstance(self.file_list_controller, LDFileListController)
        assert isinstance(self.panel, wx.Window)
        group_ctl = wx.StaticBox(panel, label = "Selection")
        sizer = wx.StaticBoxSizer(group_ctl, orient= wx.VERTICAL)
        panel.Sizer = sizer
        self.sel_action_ctrl = wx.Choice(panel, choices=[
            "Select", "Add to selection", "Remove selection"])
        sizer.Add(self.sel_action_ctrl)
        sizer.AddSpacer(2)
        sizer.Add(wx.Choice(panel, choices=["Any", "All"]))
        sub_sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(sub_sizer, 1, wx.EXPAND)
        sub_sizer.AddSpacer(10)
        sub_sizer.Add(wx.Choice(panel, choices=["File name", "Folder name", "Metadata"]))
        sub_sizer.AddSpacer(2)
        sub_sizer.Add(wx.Choice(panel, choices=["Does", "Does not"]))
        sub_sizer.AddSpacer(2)
        sub_sizer.Add(wx.Choice(panel, choices=[
            "Match regular expression", "Start with", "End with" ,"Contain"]))
        sub_sizer.AddSpacer(2)
        self.regex_ctrl = wx.TextCtrl(panel)
        sub_sizer.Add(self.regex_ctrl, 1, wx.EXPAND)
        sub_sizer.AddSpacer(2)
        sub_sizer.Add(wx.Button(panel, label="...", style=wx.BU_EXACTFIT))
        sub_sizer.AddSpacer(2)
        sub_sizer.Add(wx.Button(panel, label="+", style=wx.BU_EXACTFIT))
        sizer.AddSpacer(2)
        
        sub_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel.Sizer.Add(sub_sizer)
        self.sel_source_ctrl = wx.Choice(panel, choices=["In file list", "In DNA channel", "In GFP channel"])
        sub_sizer.Add(self.sel_source_ctrl)
        sub_sizer.AddSpacer(2)
        do_it_button = wx.Button(panel, label = "Do it")
        sub_sizer.Add(do_it_button)
        do_it_button.Bind(wx.EVT_BUTTON, self.on_action)
        for child in self.panel.GetChildren():
            if isinstance(child, wx.Choice):
                child.Select(0)
        
    def on_action(self, event):
        regex = self.regex_ctrl.Value
        selection = self.file_list_controller.file_list.get_filtered_files(
            self.filter_file)
        select_flag = self.sel_action_ctrl.StringSelection in (
            "Select", "Add to selection")
        if self.sel_action_ctrl.StringSelection == "Select":
            self.file_list_controller.select_all(False)
        self.file_list_controller.select(selection, select_flag)
        
    def filter_file(self, f):
        regex = self.regex_ctrl.Value
        return re.search(regex, f.filename) is not None
        
        
class LDFileListController(object):
    def __init__(self, file_list_ctrl, file_list):
        assert isinstance(file_list_ctrl, wx.ListCtrl)
        self.file_list_ctrl = file_list_ctrl
        assert isinstance(file_list, LDFileList)
        self.file_list = file_list
        self.item_dictionary = {}
        self.path_dictionary = {}
        self.drop_target = LDDropTarget(self.on_drop_files,
                                        self.on_drop_text)
        self.file_list_ctrl.SetDropTarget(self.drop_target)
        self.last_regexp = "^(?P<Plate>.*?)_(?P<Well>[A-Za-z]+[0-9]+)_s(?P<Site>[0-9])_w(?P<Wavelength>[0-9])\\.tif$"
        self.file_list_ctrl.Bind(wx.EVT_ERASE_BACKGROUND, self.on_erase_background)
        
    def on_erase_background(self, event):
        assert isinstance(event, wx.EraseEvent)
        dc = event.DC
        assert isinstance(dc, wx.DC)
        brush = wx.Brush(self.file_list_ctrl.GetBackgroundColour())
        dc.SetBrush(brush)
        dc.SetPen(wx.TRANSPARENT_PEN)
        width, height = self.file_list_ctrl.GetSize()
        dc.DrawRectangle(0, 0, width, height)
        if len(self.file_list) == 0:
            text = ["Drop", "files", "and","folders","here"]
            font = wx.Font(36, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL,
                           wx.FONTWEIGHT_BOLD)
            dc.SetTextForeground(wx.SystemSettings.GetColour(wx.SYS_COLOUR_GRAYTEXT))
            dc.SetFont(font)
            max_width = 0
            max_height = 0
            leading = 2
            for t in text:
                text_width, text_height = dc.GetTextExtent(t)
                max_width = max(max_width, text_width)
                max_height = max(max_height, text_height)
            total_height = max_height * len(text) + leading * (len(text)-1)
            for i, t in enumerate(text):
                dc.DrawText(t,
                            (width - max_width) / 2, 
                            (height - total_height) / 2 + max_height * i)
        
    def update_file_list(self):
        metadata_columns = list(self.file_list.get_metadata_keys())
        columns = ["File", "Channel"] + metadata_columns
        self.file_list_ctrl.ClearAll()
        self.item_dictionary = {}
        self.path_dictionary = {}
        for i, column in enumerate(columns):
            self.file_list_ctrl.InsertColumn(i, column)
        for f in self.file_list.get_files():
            assert isinstance(f, LDFile)
            fields = [f.filename, f.channel or ""]
            for i, mdkey in enumerate(metadata_columns):
                if mdkey in f.metadata:
                    fields.append(unicode(f.metadata[mdkey]))
                else:
                    fields.append("")
            idx = self.file_list_ctrl.Append(fields)
            self.item_dictionary[f.path] = idx
            self.path_dictionary[idx] = f.path

    def select(self, selection, select_state):
        count = 0
        first_item = None
        for path in sorted(selection):
            if path in self.item_dictionary:
                item = self.item_dictionary[path]
                if select_state != bool(
                    self.file_list_ctrl.IsSelected(item)):
                    if first_item is None:
                        first_item = item
                    self.file_list_ctrl.Select(item, select_state)
                    count += 1
        if select_state:
            set_status_message("Selected %d files" % count)
            if first_item is not None:
                self.file_list_ctrl.EnsureVisible(first_item)
                self.file_list_ctrl.SetFocus()
        else:
            set_status_message("Deselected %d files" % count)
                
    def get_selection(self):
        result = []
        for path, item in self.item_dictionary.iteritems():
            if self.file_list_ctrl.IsSelected(item):
                result.append(path)
        selection = LDSelection(result)
        return self.file_list.get_selection(selection)
    
    def select_all(self, select_state = True):
        self.select([f.path for f in self.file_list], select_state)
        
    def select_using_regexp(self, select_state = True):
        fn = self.get_regexp_fn()
        if fn is None:
            return
        selection = self.file_list.get_filtered_files(fn)
        self.select(selection, select_state)
        
    def get_regexp_fn(self):
        if len(self.file_list) == 0:
            wx.MessageBox("There are no files to select", caption = "No files",
                          style = wx.ICON_ERROR)
            return None
        
        regexp = edit_regexp(self.file_list_ctrl, self.last_regexp,
                             self.file_list.get_any_file().filename)
        if regexp is not None:
            self.last_regexp = regexp
            def fn(f):
                return re.search(regexp, f.filename) is not None
            return fn
        return None
    
    def add_files_by_path(self, paths):
        old_count = len(self.file_list)
        self.file_list.add_files([
            LDFile(path) for path in paths
            if path not in self.file_list.files])
        self.file_list.sort()
        self.update_file_list()
        new_count = len(self.file_list)
        set_status_message("Added %d files" % (new_count - old_count))
        
    def on_drop_files(self, x, y, filenames):
        self.add_files_by_path(filenames)
    
    def on_drop_text(self, x, y, text):
        self.on_drop_files(x, y, [t.strip() for t in text.split("\n")])
    
    def assign_channel(self, files, channel_name):
        for f in files:
            if f.path in self.item_dictionary:
                f.channel = channel_name
                idx = self.item_dictionary[f.path]
                self.file_list_ctrl.SetStringItem(idx, 1, channel_name or "")
            
class ChannelProperties(object):
    CP_MONOCHROME = "Monochrome"
    CP_COLOR = "Color"
    CP_MOVIE = "Movie"
    CP_OBJECTS = "Objects"
    CP_ILLUM = "Illumination function"
    CP_ALL = (CP_MONOCHROME, CP_COLOR, CP_MOVIE, CP_OBJECTS, CP_ILLUM)
    def __init__(self):
        self.name = "Unnamed"
        self.image_type = self.CP_MONOCHROME
        
class ChannelPropertiesDlg(wx.Dialog):
    def __init__(self, properties, *args, **kwargs):
        wx.Dialog.__init__(self, *args, **kwargs)
        self.properties = properties
        assert isinstance(properties, ChannelProperties)
        self.Sizer = wx.BoxSizer(wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.Sizer.Add(sizer, 1, wx.EXPAND)
        sizer.Add(wx.StaticText(self, label = "Name:"), 0, wx.ALIGN_LEFT)
        sizer.AddSpacer(3)
        self.name_ctrl = wx.TextCtrl(self, value = self.properties.name)
        sizer.Add(self.name_ctrl, 1, wx.EXPAND)
        self.Sizer.AddSpacer(5)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.Sizer.Add(sizer, 1, wx.EXPAND)
        sizer.Add(wx.StaticText(self, label = "Image type:"), 0, wx.ALIGN_LEFT)
        sizer.AddSpacer(3)
        self.image_type_ctrl = wx.Choice(
            self, -1, choices = properties.CP_ALL)
        self.image_type_ctrl.SetStringSelection(properties.image_type)
        sizer.Add(self.image_type_ctrl, 0, wx.ALIGN_LEFT)
        self.Sizer.AddSpacer(5)
        
        sizer = wx.StdDialogButtonSizer()
        sizer.AddButton(wx.Button(self, wx.ID_OK))
        sizer.AddButton(wx.Button(self, wx.ID_CANCEL))
        sizer.Realize()
        self.Sizer.Add(sizer, 0, wx.ALIGN_CENTER)
        
        self.Fit()
        
    def ShowModal(self):
        result = wx.Dialog.ShowModal(self)
        if result == wx.ID_OK:
            self.properties.name = self.name_ctrl.Value
            self.properties.image_type = self.image_type_ctrl.GetStringSelection()
        return result
        

class ChannelController(wx.grid.PyGridTableBase):
    def __init__(self, frame, panel, menu, file_list_controller):
        wx.grid.PyGridTableBase.__init__(self)
        self.channels = {}
        self.frame = frame
        assert isinstance(menu, wx.Menu)
        self.menu = menu
        assert isinstance(file_list_controller, LDFileListController)
        self.file_list_controller = file_list_controller
        self.panel = panel
        self.panel.Sizer = wx.BoxSizer(wx.VERTICAL)
        self.grid = wx.grid.Grid(self.panel)
        self.grid.SetTable(self)
        self.panel.Sizer.Add(self.grid, 1, wx.EXPAND)
        ch_add = wx.NewId()
        menu.Append(ch_add, "Add new channel")
        frame.Bind(wx.EVT_MENU, self.add_new_channel, id=ch_add)
        
        self.drop_target = LDDropTarget(self.on_drop_filenames,
                                        self.on_drop_text,
                                        self.on_drag_over,
                                        self.on_leave)
        self.grid.GetGridWindow().SetDropTarget(self.drop_target)
        self.panel.Sizer.Show(self.grid, False)
        self.hover_column = None
        self.panel.Bind(wx.EVT_ERASE_BACKGROUND, self.on_erase_background)
        
    def on_erase_background(self, event):
        assert isinstance(event, wx.EraseEvent)
        dc = event.DC
        assert isinstance(dc, wx.DC)
        brush = wx.Brush(self.panel.GetBackgroundColour())
        dc.SetBrush(brush)
        dc.SetPen(wx.TRANSPARENT_PEN)
        width, height = self.panel.GetSize()
        dc.DrawRectangle(0, 0, width, height)
        if not self.grid.Shown:
            text = ["Add", "Channels"]
            font = wx.Font(36, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL,
                           wx.FONTWEIGHT_BOLD)
            dc.SetTextForeground(wx.SystemSettings.GetColour(wx.SYS_COLOUR_GRAYTEXT))
            dc.SetFont(font)
            max_width = 0
            max_height = 0
            leading = 2
            for t in text:
                text_width, text_height = dc.GetTextExtent(t)
                max_width = max(max_width, text_width)
                max_height = max(max_height, text_height)
            total_height = max_height * len(text) + leading * (len(text)-1)
            for i, t in enumerate(text):
                dc.DrawText(t,
                            (width - max_width) / 2, 
                            (height - total_height) / 2 + max_height * i)
            
    def GetNumberRows(self):
        if len(self.channels) == 0:
            return 0
        return reduce(max, [len(stuff.file_list) 
                            for stuff in self.channels.values()])
    
    def GetNumberCols(self):
        return len(self.channels)
    
    def IsEmptyCell(self, row, col):
        name = self.get_channel_name_by_idx(col)
        if name is None:
            return True
        return len(self.channels[name].file_list) >= row
    
    def GetValue(self, row, col):
        name = self.get_channel_name_by_idx(col)
        stuff = self.channels[name]
        if len(stuff.file_list) <= row:
            return None
        return stuff.file_list[row].filename
    
    def GetRowLabelValue(self, row):
        return str(row+1)
        
    def GetColLabelValue(self, col):
        return self.get_channel_name_by_idx(col)
    
    def AppendCols(self, numCols):
        return True
    
    def AppendRows(self, numRows):
        return True
    
    def InsertCols(self, index, numCols):
        return True
    
    def InsertRows(self, index, numRows):
        return True
    
    def DeleteCols(self, index, numCols):
        return True
    
    def DeleteRows(self, index, numRows):
        return True
    
    def on_drop_text(self, x, y, text):
        self.on_drop_filenames(x, y, [t.strip() for t in text.split("\n")])
    
    def on_drop_filenames(self, x, y, filenames):
        x, y = self.grid.CalcUnscrolledPosition(x, y)
        col = self.grid.XToCol(x)
        name = self.get_channel_name_by_idx(col)
        if name is not None:
            self.file_list_controller.add_files_by_path(filenames)
            selection = LDSelection(filenames)
            files = self.file_list_controller.file_list.get_selection(selection)
            self.add_files_to_channel(name, files)
            self.grid.DeselectCol(col)
    
    def on_drag_over(self, x, y, d):
        x, y = self.grid.CalcUnscrolledPosition(x, y)
        old_hover = self.hover_column
        self.hover_column = self.grid.XToCol(x)
        if self.hover_column >= 0 and self.hover_column < len(self.channels):
            return wx.DragCopy
        return wx.DragNone
    
    def on_leave(self):
        self.hover_column = None
        
    class ChannelStuff(object):
        def __init__(self, name, properties, file_list, column_idx, menu_item):
            self.name = name
            self.properties = properties
            self.file_list = file_list
            self.column_idx = column_idx
            self.menu_item = menu_item
            
    def get_channel_name_by_idx(self, idx):
        channel_names = filter(
            (lambda name: self.channels[name].column_idx == idx),
            self.channels.keys())
        return None if len(channel_names) == 0 else channel_names[0]
            
    def add_new_channel(self, event):
        if not self.grid.Shown:
            self.panel.Sizer.Show(self.grid, True)
            self.panel.Layout()
        new_channel_properties = ChannelProperties()
        with ChannelPropertiesDlg(new_channel_properties, self.frame) as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                column = self.grid.GetNumberCols()
                name = new_channel_properties.name
                add_selection_id = wx.NewId()
                menu = wx.Menu()
                menu.Append(add_selection_id, "Add files from file list")
                def add_fn(event):
                    self.add_files_from_file_list(new_channel_properties.name)
                self.frame.Bind(wx.EVT_MENU, add_fn, id = add_selection_id)
                del_selection_id = wx.NewId()
                menu.Append(del_selection_id, "Delete channel")
                def del_fn(event):
                    self.delete_channel(new_channel_properties.name)
                self.frame.Bind(wx.EVT_MENU, del_fn, id = del_selection_id)
                remove_selection_id = wx.NewId()
                def remove_fn(event):
                    self.remove_selection(new_channel_properties.name)
                menu.Append(remove_selection_id, "Remove selected files")
                self.frame.Bind(wx.EVT_MENU, remove_fn, id=remove_selection_id)
                properties_id = wx.NewId()
                menu.Append(properties_id, "Properties")
                def props_fn(event):
                    self.edit_properties(new_channel_properties.name)
                self.frame.Bind(wx.EVT_MENU, props_fn, id = properties_id)
                menu_item = self.menu.AppendSubMenu(menu, name)
                self.channels[name] = self.ChannelStuff(
                    name, new_channel_properties, LDFileList(), column, menu_item)
                self.grid.AppendCols(1)
                self.grid.SetColLabelValue(column, name)
                tm = wx.grid.GridTableMessage(
                    self,
                    wx.grid.GRIDTABLE_NOTIFY_COLS_INSERTED,
                    0, 1)
                self.grid.ProcessTableMessage(tm)
                self.grid.ForceRefresh()
                set_status_message('Added channel, "%s"' % name)
    
    def add_files_from_file_list(self, channel_name):
        new_files = self.file_list_controller.get_selection()
        self.add_files_to_channel(channel_name, new_files)
        
    def add_files_to_channel(self, channel_name, new_files):
        channel_stuff = self.channels[channel_name]
        file_list = channel_stuff.file_list
        old_len = self.GetNumberRows()
        assert isinstance(file_list, LDFileList)
        for f in new_files:
            if f.channel is not None and f.channel != channel_name:
                self.channels[f.channel].file_list.remove_files((f,))
                
        file_list.add_files(new_files)
        new_len = self.GetNumberRows()
        if old_len < new_len:
            tm = wx.grid.GridTableMessage(
                self,
                wx.grid.GRIDTABLE_NOTIFY_ROWS_INSERTED,
                old_len, new_len-old_len)
            self.grid.ProcessTableMessage(tm)
        self.grid.ForceRefresh()
        self.file_list_controller.assign_channel(new_files, channel_name)
        set_status_message("Added %d files to %s channel" % (len(new_files),
                                                             channel_name))
            
    def delete_channel(self, channel_name):
        channel_stuff = self.channels[channel_name]
        self.file_list_controller.assign_channel(channel_stuff.file_list, None)
        idx = channel_stuff.column_idx
        for other_stuff in self.channels.values():
            if other_stuff.column_idx > idx:
                other_stuff.column_idx -= 1
        del self.channels[channel_name]
        self.menu.RemoveItem(channel_stuff.menu_item)
        set_status_message("Channel %s deleted" % channel_name)
        tm = wx.grid.GridTableMessage(
            self,
            wx.grid.GRIDTABLE_NOTIFY_COLS_DELETED,
            idx, 1)
        self.grid.ProcessTableMessage(tm)
        self.grid.ForceRefresh()
        if len(self.channels) == 0:
            self.panel.Sizer.Show(self.grid, False)
            self.panel.Layout()
            
    def remove_selection(self, channel_name):
        old_len = self.GetNumberRows()
        rows = self.grid.GetSelectedCells()
        channel_stuff = self.channels[channel_name]
        del channel_stuff.file_list[rows]
        new_len = self.GetNumberRows()
        if old_len > new_len:
            tm = wx.grid.GridTableMessage(
                self,
                wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED,
                new_len, old_len-new_len)
            self.grid.ProcessTableMessage(tm)
            
        self.grid.ForceRefresh()
            
    
    def edit_properties(self, channel_name):
        channel_stuff = self.channels[channel_name]
        channel_properties = channel_stuff.properties
        with ChannelPropertiesDlg(channel_properties, self.frame) as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                if channel_name != channel_properties.name:
                    self.file_list_controller.assign_channel(
                        channel_stuff.file_list, channel_properties.name)
                    self.menu.SetLabel(channel_stuff.menu_item.Id,
                                       channel_properties.name)
                    del self.channels[channel_name]
                    self.channels[channel_properties.name] = channel_stuff
                    self.grid.ForceRefresh()
                set_status_message("Channel %s properties edited" % 
                                   channel_properties.name)
status_bar = None                
def main():
    global status_bar
    app = wx.PySimpleApp(True)
    frame = wx.Frame(None, size=(800,600))
    sizer = wx.BoxSizer(wx.VERTICAL)
    frame.Sizer = sizer
    select_panel = wx.Panel(frame)
    sizer.Add(select_panel, 0, wx.EXPAND)
    splitter = wx.SplitterWindow(frame)
    sizer.Add(splitter, 1, wx.EXPAND)
    
    file_list_panel = wx.Panel(splitter)
    file_list_panel.Sizer = wx.BoxSizer(wx.VERTICAL)
    file_list_panel.Sizer.Add(wx.StaticText(file_list_panel,
                                            label = "File List"),
                              0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5)
    file_list_ctrl = wx.ListCtrl(file_list_panel,
                                 style = wx.LC_REPORT)
    file_list_panel.Sizer.Add(file_list_ctrl, 1, wx.EXPAND)
    file_list = LDFileList()
    file_list_controller = LDFileListController(file_list_ctrl, file_list)
    
    right_panel = wx.Panel(splitter)
    right_panel.Sizer = wx.BoxSizer(wx.VERTICAL)
    channel_panel = wx.Panel(splitter)
    status_bar = wx.StatusBar(frame)
    status_bar.SetFieldsCount(1)
    frame.Sizer.Add(status_bar, 0, wx.EXPAND)
    
    menu_bar = wx.MenuBar()
    file_menu = wx.Menu()
    file_menu.Append(wx.ID_OPEN, "Open")
    file_menu.Append(wx.ID_SAVE, "Save")
    menu_bar.Append(file_menu, "File")
    edit_menu = wx.Menu()
    select_all_id = wx.NewId()
    frame.Bind(wx.EVT_MENU, lambda event: file_list_controller.select_all(True),
               id = select_all_id)
    edit_menu.Append(select_all_id, "Select all")
    
    deselect_all_id = wx.NewId()
    frame.Bind(wx.EVT_MENU, lambda event: file_list_controller.select_all(False),
               id = deselect_all_id)
    edit_menu.Append(deselect_all_id, "Deselect all")
    fl_select_regexp_id = wx.NewId()
    frame.Bind(
        wx.EVT_MENU, 
        lambda event: file_list_controller.select_using_regexp(True),
        id = fl_select_regexp_id)
    fl_deselect_regexp_id = wx.NewId()
    frame.Bind(
        wx.EVT_MENU, 
        lambda event: file_list_controller.select_using_regexp(False),
        id = fl_deselect_regexp_id)
    re_select_menu = wx.Menu()
    re_select_menu.Append(fl_select_regexp_id, "In file list")
    edit_menu.AppendMenu(-1,"Regexp select", re_select_menu)
    re_deselect_menu = wx.Menu()
    re_deselect_menu.Append(fl_deselect_regexp_id, "In file list")
    edit_menu.AppendMenu(-1,"Regexp deselect", re_deselect_menu)
    menu_bar.Append(edit_menu, "Edit")
    
    channel_menu = wx.Menu()
    channel_controller = ChannelController(frame, channel_panel, channel_menu,
                                           file_list_controller)
    select_controller = LDSelectionUIController(select_panel, file_list_controller)
    menu_bar.Append(channel_menu, "Channels")
    frame.SetMenuBar(menu_bar)

    splitter.SplitVertically(file_list_panel, channel_panel)
    
    frame.Layout()
    frame.Show()
    
    app.MainLoop()
    
def set_status_message(text):
    status_bar.SetStatusText(text, 0)
    
if __name__ == "__main__":
    main()