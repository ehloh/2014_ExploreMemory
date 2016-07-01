function f_sendemail(id,subject,message)
% Original Smittenaar: f_sendmail(id,subject,message,attachment)
%% SEND_MAIL_MESSAGE send email to gmail once calculation is done
% Example
% send_mail_message('petersmittenaar','this is the subject','This is the main message','results.mat')
% will send email to petersmittenaar@gmail.com, with results.doc attached.

% Pradyumna
% June 2008
% adapted by Peter Smittenaar, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Your gmail ID and password
%(from which email ID you would like to send the mail)
switch id
    case 'learnreward'
        mail = 'learnreward@gmail.com';    %Your GMail email address
        password = 'memexpts';          %Your GMail password
    case 'kurzlich'
        mail='kurzlich@gmail.com';
        password='koko9090';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    id = 'learnreward';
    message = 'Subject is done';
    subject = 'subject';
    attachment = [];
elseif nargin == 1
    message = 'Go and look at your results';
    subject = subject;
    attachment = [];
elseif nargin == 2
    message = 'Go and look at your results';
    attachment = [];
elseif nargin == 3
    attachment = [];
end

% Send Mail ID
emailto = strcat(id,'@gmail.com');
%% Set up Gmail SMTP service.
% Then this code will set up the preferences properly:
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);

% Gmail server.
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%% Send the email
if strcmp(mail,'GmailId@gmail.com')
    disp('Please provide your own gmail.')
    disp('You can do that by modifying the first two lines of the code')
    disp('after the comments.')
end

if isempty(attachment)
    sendmail(emailto,subject,message)
else
    sendmail(emailto,subject,message)
end
end
